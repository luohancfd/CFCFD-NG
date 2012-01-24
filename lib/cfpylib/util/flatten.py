# flatten.py
#
# Taken from the Python Cookbook
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/363051
# PJ, 03-May-2007
#
def flatten(sequence):
    def rflat(seq2):
        seq = []
        for entry in seq2:
            if contains_sequence([entry]):
                seq.extend([i for i in entry])
            else:
                seq.append(entry)
        return seq

    def contains_sequence(sequence):
        for i in sequence:
            ## all sequences have '__contains__' in their dir()
            if ('__contains__' in dir(i) and type(i) != str and type(i) != dict):
                return True
        return False

    seq = [sequence][:]    ## in case parameter isn't already a sequence
    while contains_sequence(seq):
        seq = rflat(seq)
    return seq


if __name__ == '__main__':
    print "Try flattening sequences..."
    a = [1,2,[3,4,[5,]]]
    print "a=", a, "flatten(a)=", flatten(a)
    a = [(3,4,[5,]),6]
    print "a=", a, "flatten(a)=", flatten(a)
    a = [(3,4,[5,]),[]]
    print "a=", a, "flatten(a)=", flatten(a)
    a = [1,2,]
    print "a=", a, "flatten(a)=", flatten(a)
    a = []
    print "a=", a, "flatten(a)=", flatten(a)
    print "Done"
