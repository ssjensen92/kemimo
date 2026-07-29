from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ModelConfig:
    datadir: str
    nlayers: int
    layer_thickness: float
    ignore_missing_h: bool
    ignore_missing_ebind: bool
    all_species_to_dust: bool
    limit_2body: bool
    respect_gasphase_limits: bool
    multiprocessing: bool
    do_swap: bool
    h2_spin: bool
    encounter_desorption: bool
    reaction_diffusion_competition: bool
    thin_ice_approximation: bool


_FIELD_TYPES = {
    "datadir": str,
    "nlayers": int,
    "layer_thickness": float,
    "ignore_missing_h": bool,
    "ignore_missing_ebind": bool,
    "all_species_to_dust": bool,
    "limit_2body": bool,
    "respect_gasphase_limits": bool,
    "multiprocessing": bool,
    "do_swap": bool,
    "h2_spin": bool,
    "encounter_desorption": bool,
    "reaction_diffusion_competition": bool,
    "thin_ice_approximation": bool,
}


def _strip_comment(line):
    quote = None
    result = []
    for character in line:
        if character in ("'", '"'):
            if quote is None:
                quote = character
            elif quote == character:
                quote = None
        if character == "!" and quote is None:
            break
        result.append(character)
    return "".join(result)


def _split_assignments(body):
    quote = None
    item = []
    for character in body:
        if character in ("'", '"'):
            if quote is None:
                quote = character
            elif quote == character:
                quote = None
        if character in (",", "\n") and quote is None:
            assignment = "".join(item).strip()
            if assignment:
                yield assignment
            item = []
        else:
            item.append(character)
    assignment = "".join(item).strip()
    if assignment:
        yield assignment


def _parse_bool(name, value):
    normalized = value.strip().lower()
    if normalized in ("1", ".true.", "true"):
        return True
    if normalized in ("0", ".false.", "false"):
        return False
    raise ValueError(
        "%s must be 1 or 0 (Fortran .true./.false. are also accepted)" % name)


def _parse_value(name, value, expected_type):
    value = value.strip()
    if expected_type is bool:
        return _parse_bool(name, value)
    if expected_type is str:
        if len(value) < 2 or value[0] not in ("'", '"') or value[-1] != value[0]:
            raise ValueError("%s must be a quoted string" % name)
        return value[1:-1]
    try:
        return expected_type(value.replace("d", "e").replace("D", "e"))
    except ValueError as exc:
        raise ValueError(
            "invalid value for %s: %s" % (name, value)) from exc


def load_model_config(filename="config.nml"):
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError("model configuration not found: %s" % path)

    lines = [_strip_comment(line) for line in path.read_text().splitlines()]
    content = "\n".join(lines).strip()
    if not content.lower().startswith("&kemimo"):
        raise ValueError("%s must start with the &kemimo namelist" % path)

    body = content[len("&kemimo"):].strip()
    if not body.endswith("/"):
        raise ValueError("%s must end the &kemimo namelist with /" % path)
    body = body[:-1]

    values = {}
    for assignment in _split_assignments(body):
        if "=" not in assignment:
            raise ValueError(
                "invalid namelist assignment in %s: %s" % (path, assignment))
        name, value = assignment.split("=", 1)
        name = name.strip().lower()
        if name not in _FIELD_TYPES:
            raise ValueError(
                "unknown model configuration variable: %s" % name)
        if name in values:
            raise ValueError(
                "duplicate model configuration variable: %s" % name)
        values[name] = _parse_value(name, value, _FIELD_TYPES[name])

    missing = [name for name in _FIELD_TYPES if name not in values]
    if missing:
        raise ValueError(
            "missing model configuration variable(s): %s"
            % ", ".join(missing))

    config = ModelConfig(**values)
    if config.nlayers not in (0, 1, 2):
        raise ValueError("nlayers must be 0, 1, or 2")
    if config.layer_thickness <= 0.0:
        raise ValueError("layer_thickness must be positive")
    return config
