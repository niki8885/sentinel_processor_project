STAC_API_URL = "https://earth-search.aws.element84.com/v1"

DEFAULT_BBOX_HALF_DEG: float = 0.05
DEFAULT_MAX_ITEMS: int   = 50
DEFAULT_KEEP_ITEMS: int   = 10
DEFAULT_LOOKBACK_DAYS: int   = 30

DEFAULT_OUTPUT_DIR: str  = "data"
DEFAULT_NAME_TEMPLATE: str = "{date}_{lat}_{lon}"
SAVE_VALIDATION_REPORT: bool = True

MAX_CLOUD_THRESHOLD: float = 0.30
MIN_CONFIDENCE_TO_SAVE: float = 0.01