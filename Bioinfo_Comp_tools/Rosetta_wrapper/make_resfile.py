def resfile(mutations:list, filename:str, workdir:str):
    """
    generate resfile for rosetta flex_ddg given modern style input
    =======================================
    mutation here has diferent definition, the format is <chain><index><mutation resname>
    """
    with open(f"{workdir}/{filename}", 'w') as f:
        f.write("NATAA\nstart\n")

        for mut in mutations:
            f.write(f"{mut[1:-1]} {mut[0]} PIKAA {mut[-1]}\n")

    return f"{workdir}/{filename}"