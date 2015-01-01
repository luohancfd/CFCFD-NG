def set_perturbed_values(NominalValue, PerturbationMagnitude, TypeOfPerturbation, Levels):
    #
    if TypeOfPerturbation == "relative":
        delta = [k/100.0*NominalValue for k in PerturbationMagnitude]
    else:
        delta = PerturbationMagnitude
    #
    if Levels == 3:
        if len(delta) == 1:
            # have [delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[0]]
        elif len(delta) == 2:
            # have [+delta, -delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[1]]
    elif levels == 5:
        if len(delta) == 1:
            # have [delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[0],\
                                    NominalValue+2.0*delta[0],\
                                    NominalValue-2.0*delta[0]]
        elif len(delta) == 2:
            # have [delta, "2"*delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[0],\
                                    NominalValue+delta[1],\
                                    NominalValue-delta[1]]
        elif len(delta) == 3:
            # have [+delta, -delta, "2"*delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[1],\
                                    NominalValue+delta[2],\
                                    NominalValue-delta[2]]
        elif len(data) == 5:
            # have [+delta, -delta, +"2"*delta, -"2"*delta]
            values = [NominalValue, NominalValue+delta[0],\
                                    NominalValue-delta[1],\
                                    NominalValue+delta[2],\
                                    NominalValue-delta[3]]
        
    return values
