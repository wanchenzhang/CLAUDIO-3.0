"""Keep every M4 candidate alongside its selection history."""
import json


def build_candidate_audit(candidates, selected, selection_log):
    """Describe final-table selection without changing candidate source fields."""
    audit_columns = (
        "m4_selection_status", "m4_exclusion_events", "m4_audit_rule_version",
    )
    for column in audit_columns:
        if column in candidates.columns:
            raise ValueError(f"Reserved audit column already exists: {column}")

    labels = candidates.index.map(str)
    selected_labels = selected.index.map(str)
    if not labels.is_unique or not selected_labels.is_unique:
        raise ValueError("Candidate audit indices must be unique as text")

    all_indices = set(labels)
    retained = set(selected_labels)
    if not retained.issubset(all_indices):
        raise ValueError("Unknown selected candidate in M4 audit")

    events_by_index = {index: [] for index in labels}
    excluded = all_indices - retained
    for event in selection_log:
        index = str(event["removed_index"])
        if index not in excluded:
            raise ValueError(f"Invalid exclusion event for candidate {index}")
        if any(str(value) not in all_indices
               for value in event["preferred_indices"]):
            raise ValueError(f"Unknown preferred candidate for {index}")
        events_by_index[index].append(event)

    logged = {index for index, events in events_by_index.items() if events}
    if logged != excluded:
        raise ValueError("An excluded candidate has no exclusion event")

    result = candidates.copy(deep=True)
    result["m4_selection_status"] = [
        "retained" if index in retained else "excluded" for index in labels
    ]
    result["m4_exclusion_events"] = [
        json.dumps(events_by_index[index], ensure_ascii=False, allow_nan=False)
        for index in labels
    ]
    result["m4_audit_rule_version"] = "m4-candidate-audit-v1"
    return result