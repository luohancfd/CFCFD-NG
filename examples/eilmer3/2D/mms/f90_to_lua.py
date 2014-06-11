#! /usr/bin/env python
# Translate f90 code for evaluating source-term expressions
# into Lua-compatible code, then paste it into the Lua file
# that gets called by Eilmer3.
# PJ, 29-May-2011

fin = open('source_terms.f90', 'r')
text = fin.read()
fin.close()
text = text.replace('&\n', '')
text = text.replace('&', '')
text = text.replace('%pi', 'math.pi')
text = text.replace('sin', 'math.sin')
text = text.replace('cos', 'math.cos')
text = text.replace('exp', 'math.exp')
source_terms_text = text.replace('**', '^')

fin = open('udf-source-template.lua', 'r')
template_text = fin.read()
fin.close()
lua_text = template_text.replace('<insert-source-terms-here>',
                                 source_terms_text)

fout = open('udf-source.lua', 'w')
fout.write(lua_text)
fout.close()
print 'Done converting to Lua.'


