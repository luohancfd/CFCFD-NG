
module(..., package.seeall)

-----
-- Adapted from Listing 12.2 in
-- Programming in Lua, 2nd edition
-----
function serialise(o, f, indent)
   indent = indent or ""
   if type(o) == "number" then
      f:write(o)
   elseif type(o) == "string" then
      f:write(string.format("%q", o))
   elseif type(o) == "table" then
      local old_indent = indent
      indent = indent.."  "
      f:write("{\n")
      for k,v in pairs(o) do
	 if type(k) == "string" then
	    f:write(string.format("%s%s = ", indent, k))
	 else
	    f:write(string.format("%s", indent))
	 end
	 serialise(v, f, indent)
	 f:write(",\n")
      end
      indent = old_indent
      f:write(string.format("%s}", indent))
   else
      error("cannot serialise a " .. type(o))
   end
end
