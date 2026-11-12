"""Allow the past-only boundary of PP without changing the archived S3 engine."""
from dataclasses import dataclass, asdict
from model_s3 import Parameters as BaseParameters, simulate, weights, seed_for


@dataclass(frozen=True)
class Parameters(BaseParameters):
    def __post_init__(self):
        if self.kernel == 'PP' and self.a == 0:
            # Check every other parameter through the unchanged S3 validator.
            BaseParameters(**(asdict(self) | {'a': 0.5}))
        else:
            super().__post_init__()
