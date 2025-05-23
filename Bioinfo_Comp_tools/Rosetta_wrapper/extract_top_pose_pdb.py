def extract(score_file, topk):
    with open(score_file, "r") as fh:
        scores = fh.readlines()[2:]
    return