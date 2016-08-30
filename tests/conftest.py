import pytest

@pytest.fixture
def config_path():
    """ Return the path to the config yaml file """
    return "tests/fixtures/config.yaml"

@pytest.yield_fixture(scope='function')
def config_stream(config_path):
    """ yield a stream to the config file """
    stream = open(config_path, 'r')
    yield stream
    stream.close()
