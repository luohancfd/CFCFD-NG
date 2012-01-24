from spark import GenericScanner, GenericParser

class Token:
    def __init__(self, type=type, attr=None):
	self.type = type
	self.attr = attr

    def __cmp__(self, o):
	return cmp(self.type, o)

    def __repr__(self):
	if self.attr == None:
	    return self.type
	else:
	    return self.type+"("+self.attr+")"

class _ReactionScanner(GenericScanner):
    
    def __init__(self):
        GenericScanner.__init__(self)
    
    def tokenize(self, input):
        self.rv = []
        GenericScanner.tokenize(self, input)
        return self.rv
    
    def t_whitespace(self, s):
        r' \s+ '
        pass

    def t_number(self, s):
        r' \d+ '
        t = Token(type='number', attr=int(s))
        self.rv.append(t)

    def t_plus(self, s):
        r' \+ '
        self.rv.append(Token(type=s))

    def t_arrows(self, s):
        r' <=> | =>'
        self.rv.append(Token(type=s))

    def t_bracket(self, s):
        r' \( | \) '
        self.rv.append(Token(type=s))

class ReactionScanner(_ReactionScanner):
    def __init__(self):
        _ReactionScanner.__init__(self)

    def t_chemical(self, s):
        r' ([A-Z][a-z]?[0-9]*)+[+\-]? | e\- '
        self.rv.append(Token(type='chemical', attr=s))

class ReactionParser(GenericParser):
    def __init__(self, start='mechanism'):
        GenericParser.__init__(self, start)

    def p_mechanism_1(self, args):
        ' mechanism ::= reaction <=> reaction '
        return {'reactants': args[0],
                'products': args[2],
                'direction': "both_ways"}

    def p_mechanism_2(self, args):
        ' mechanism ::= reaction => reaction '
        return {'reactants': args[0],
                'products': args[2],
                'direction': "forward_only"}

    def p_reaction_1(self, args):
        ' reaction ::= reaction + participant '
        args[0].append(args[2])
        return args[0]

    def p_reaction_2(self, args):
        ' reaction ::= reaction ( + participant )'
        args[3]['pd_term'] = True
        args[0].append(args[3])
        return args[0]

    def p_reaction_3(self, args):
        ' reaction ::= participant '
        return [args[0]]

    def p_participant_1(self, args):
        ' participant ::= chemical '
        return {'coeff': 1, 'sp': args[0].attr}

    def p_participant_2(self, args):
        ' participant ::= number chemical '
        return {'coeff': args[0].attr, 'sp': args[1].attr}

def parse_reaction(equation):
    scanner = ReactionScanner()
    parser = ReactionParser()
    return parser.parse(scanner.tokenize(equation))

if __name__ == '__main__':
    print "Testing parse_reaction..."
    input = "CH3 + H ( + M ) <=> CH2 + 2 H (+ M)"
    scanner = ReactionScanner()
    tokens = scanner.tokenize(input)
    print "Scanned tokens:"
    print tokens
    parser = ReactionParser()
    result = parser.parse(tokens)
    print "Grammar parsing:"
    print result
