SIMSI_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(SIMSI_DIR)MakefileShared

dependencies:
	git config --global credential.helper cache

jump:
	$(DOCKER_CMD) \
		$(IMAGE) bash

# --no-cache
build: dependencies
	docker build -f Dockerfile -t $(IMAGE) . || (exit 1)

package_build:
	poetry build

publish_test:
	poetry config repositories.testpypi https://test.pypi.org/legacy/
	poetry publish -r testpypi

publish:
	poetry publish


run_simsi_transfer: rm_err_file
	$(DOCKER_CMD) \
		$(IMAGE) python3 -u -m simsi_transfer \
			--mq_txt_folder /root/data/combined/txt \
			--raw_folder /root/data \
			--output_folder /root/data/simsi_output \
			--stringencies 20,15 \
			--tmt_requantify \
			${EXTRA_ARGS} \
			|| (echo "2" > $(DATA)err.out; exit 2)

compress: run_simsi_transfer
	zip -j -r -9 "$(DATA)/results.zip" "$(DATA)/simsi_output/summaries" || (echo "3" > $(DATA)err.out; exit 3)

all: compress


run_local:
	python3 -u -m simsi_transfer "$(DATA)"

clean_data_folder: 
	rm -rf "$(DATA)/{proc,msms,percolator,mzML,msms.prosit}"


integration_test:
	python3 -um simsi_transfer \
		--mq_txt_folder ${DATA}/combined/txt/ \
		--raw_folder ${DATA} \
		--output_folder ${DATA}/simsi_output \
		--num_threads 3 \
		--stringencies 20,15 \
		${EXTRA_ARGS}

integration_test_meta:
	python3 -um simsi_transfer \
		--meta_input_file ${DATA}/meta_input_file.txt \
		--output_folder ${DATA}/simsi_output \
		--num_threads 3 \
		--stringencies 20,15 \
		--tmt_requantify \
		${EXTRA_ARGS}

test:
	python3 -m pytest --cov=simsi_transfer --cov-report html --cov-report term tests/unit_tests
