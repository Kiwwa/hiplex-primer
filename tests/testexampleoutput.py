import difflib

print '-' * 80
print "Executing compare between generated output log and expected output log..."
print '-' * 80

example_output = open('../exampledata/output/output.log', 'r').read()
static_compare = open('file_compare/static_compare_output.log', 'r').read()

diff2 = difflib.context_diff(example_output.splitlines(), static_compare.splitlines())

if list(diff2) == []:
	print "No differences between 'generated example output log' and 'expected example output log'."
else:
	print '\n'.join(list(diff2))