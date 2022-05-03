echo "GitHub release downloads"
curl -s -i https://api.github.com/repos/kusterlab/simsi-transfer/releases -H "Accept: application/vnd.github.manifold-preview+json" | grep "\"name\"\|\"download_count\""
#https://pypistats.org/api/packages/simsi_transfer/recent

echo "PyPI downloads"
curl -s -i https://pypistats.org/api/packages/simsi_transfer/recent -H "Accept: application/vnd.github.manifold-preview+json" | grep "\"last_day\""
