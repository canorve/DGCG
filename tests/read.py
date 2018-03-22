

    with open(File) as INCONF:

        # All lines including the blank ones
        lines = (line.rstrip() for line in INCONF)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        for line in lines:

            (params) = line.split()

            if params[0] == "Img":     # input image
                try:
                    (param, Img) = line.split()
                    parvar.Img = str(Img)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SexCat":     # Sextractor catalog
                try:
                    (param, SexCat) = line.split()
                    parvar.SexCat = str(SexCat)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

