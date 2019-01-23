import os


def get_project_root():
    path = os.path.dirname(__file__)

    return os.path.abspath(os.path.join(path, os.pardir))


if __name__ == "__main__":
    print(get_project_root())
