#!/usr/bin/env/python

import unittest
import primer_design

def isodd(n):
	return n % 2 == 1

class isoddTests(unittest.TestCase):

	def testone(self):
		self.failUnless(isodd(1))

	def testtwo(self):
		self.failIf(isodd(1))


def main():
	unittest.main()

if __name__ == '__main__':
	main()