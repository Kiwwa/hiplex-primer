import difflib

print '\n', '-' * 80
print "Executing compare between generated output log and expected output log..."
print '-' * 80

example_output = open('../exampledata/exampleoutput/output.log', 'r').read()
static_compare = open('file_compare/static_compare_output.log', 'r').read()

diff2 = list(difflib.context_diff(example_output.splitlines(), static_compare.splitlines()))

if diff2 == []:
	print "No differences between current output and expected output; this is good.", "\n"
else:
	print '\n'.join(list(diff2)), "\n"
