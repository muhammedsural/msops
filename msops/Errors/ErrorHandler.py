class KeyError(Exception):
    """Custom error that is raised when key is not available."""

    def __init__(self, list_keys: list[str], target_keys: str, message: str) -> None:
        self.requested_days = list_keys
        self.remaining_days = target_keys
        self.message = message
        super().__init__(message)

class ShapeError(Exception):
    """Custom error that is raised when not enough vacation days are available."""
    def __init__(self, future_index: int, target_index: int, message: str) -> None:
        self.first_shape = future_index
        self.second_shape = target_index
        self.message = message
        super().__init__(message)

