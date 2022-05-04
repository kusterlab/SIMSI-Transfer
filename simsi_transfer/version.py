import re
import os


def get_version_from_pyproject():
    """
    Implements a getting version from file flow for developers,
    without installed package
    """
    compiled_version_regex = re.compile(r"\s*version\s*=\s*[\"']\s*([-.\w]{3,})\s*[\"']\s*")
    file_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/pyproject.toml"
    try:
        with open(file_path, mode="r", encoding="utf-8") as file_with_version:
            for line in file_with_version:
                version = compiled_version_regex.search(line)  # type: ignore # noqa: E501
                if version is not None:
                    return version.group(1).strip()
    except (FileNotFoundError, UnicodeDecodeError):
        pass
    return None
