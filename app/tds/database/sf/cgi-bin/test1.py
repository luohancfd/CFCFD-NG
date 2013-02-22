#!/usr/bin/env python

# import cgi

text_to_send = """<html>
<head><title>test 1</title></head>

<body>
<H2>Something big again</H2>
something else...
</body>
</html>
"""

print "Content-type: text/html"
print "Content-length: ", len(text_to_send)
print ""
print text_to_send,
