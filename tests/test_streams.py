from bionumpy.streams.left_join import left_join
import pytest 


@pytest.fixture
def left_stream():
    return zip(['0', '1', '2'], [0, 1, 2])


@pytest.fixture
def right_stream():
    return zip(['0', '2'], [10, 20])


@pytest.fixture
def full_right_stream():
    return zip(['0', '1','2'], [10, 15, 20])


@pytest.fixture
def leftover_right_stream():
    return zip(['0', '1','2', '3'], [10, 15, 20, 30])


def test_left_join(left_stream, right_stream):
    assert list(left_join(left_stream, right_stream)) == [
        ('0', 0, 10), ('1', 1, None), ('2', 2, 20)]


def test_left_join_full(left_stream, full_right_stream):
    assert list(left_join(left_stream, full_right_stream)) == [
        ('0', 0, 10), ('1', 1, 15), ('2', 2, 20)]


def test_left_join_leftover(left_stream, leftover_right_stream):
    with pytest.raises(Exception):
        list(left_join(left_stream, leftover_right_stream))
