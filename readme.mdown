random-js
=========

A fairly close port of the Python Standard Library’s random module ([docs][pyranddocs], [source][pyrandsource]), but using a fairly simple [multiply with carry][mwc] + [xorshift][] PRNG, instead of the Mersenne Twister that Python uses by default. This PRNG should be suitable for most Monte Carlo simulations likely to run in a browser, or for purposes like procedural art. A slower but higher quality PRNG is also provided in the HighQualityRandom class. Only a few methods need be overridden to add a custom PRNG, if this still doesn’t do the trick.

  [pyranddocs]: http://docs.python.org/py3k/library/random.html
  [pyrandsource]: http://hg.python.org/cpython/file/tip/Lib/random.py
  [mwc]: http://en.wikipedia.org/wiki/Multiply-with-carry
  [xorshift]: http://en.wikipedia.org/wiki/Xorshift

Like the Python version, this module provides useful functions for dealing both with choosing integers/array items, and for generating random numbers following several common statistical distributions.

```coffeescript
random = new require('random-py').Random

# set the PRNG’s seed, so that the precise sequence of numbers can be
# regenerated. When the class is first constructed, the PRNG is seeded
# with `+new Date`
random.seed(12345)

# choose a floating point number in the range [0, 1)
random.random()

# choose a floating point number in the range [1.5, 10)
random.uniform(1.5, 10)

# choose an integer N in the range 2 <= N < 50, by 2
random.randrange(2, 50, 2)

# or choose an integer N in the range 0 <= N < 45
random.randrange(45)

# choose an integer N in the range 5 <= N <= 18, endpoint included
random.randint(5, 18)

# randomly choose an element of an array
array = 'abcdefg'.split('')
random.choice(array)

# choose 4 elements from the array, ordered, chosen without replacement
random.sample(array, 4)

# randomly shuffle the array, in place
random.shuffle(array)

# choose a random number from the standard normal distribution
random.gauss()

# choose a random number from the normal distribution with mean 5 and
# standard deviation 5
random.gauss(5, 5)

# choose from the triangular distribution on range [10, 20) with
# mode (peak) 18
random.triangular(10, 20, 18)

# choose from the triangular distribution on range [0, 1) with mode 0.5
random.triangular()

# choose from the log normal distribution. the log of this distribution
# is the normal distibution with mean 5 and standard deviation 5
random.lognormvariate(5, 5)

# choose from the Von Mises distribution, an analog of the normal distribution
# wrapped around a circle, with mean angle π, and concentration parameter π/2
random.vonmisesvariate(Math.PI, Math.PI/2)

# other distributions:
#   - expovariate
#   - gammavariate
#   - betavariate
#   - paretovariate
#   - weibullvariate

# it’s possible to save and restore the state of the PRNG, allowing the same
# set of random numbers to be generated in the same order:
some_state = random.getstate()

a = random.random()
random.random() for i in [0...1000]

random.setstate(some_state)

random.random() == a

# * * * * * * * * * * * * *

# If the built-in PRNG doesn’t meet your needs, it is easy to
# override with your own PRNG. But this module also ships with
# a couple of alternatives.
{BaseRandom, BuiltinRandom, HighQualityRandom} = require('random-py')

# First, `BuiltinRandom` generates random numbers approximately 10
# times as fast as `Random`. It calls the built-in `Math.random`
# function twice to construct a random number between 0 and 1
# with a full 52 bits of entropy. Unfortunately, typical JavaScript
# engines have PRNGs with rather poor performance on statistical
# tests of randomness, and this class also does not support setting
# a custom seed or saving/restoring the PRNG state.
random = new BuiltinRandom
random.random()

# The other PRNG provided has a much longer period and should pass
# more rigorous statistical tests, at the cost of running roughly 8–10
# times slower:
random = new HighQualityRandom
random.random()

# It is also quite straight-forward to implement a better custom PRNG:
class XKCDRandom extends BaseRandom
    # http://xkcd.com/221/
    _randint32: -> 4
    _seed: (j) -> # ignore j
    _getstate: ->
    _setstate: ->

class DilbertRandom extends BaseRandom
    # http://dilbert.com/fast/2001-10-25/
    _randint32: -> 9
    _seed: (j) -> # ignore j
    _getstate: ->
    _setstate: ->
```

[![Random Number](http://imgs.xkcd.com/comics/random_number.png)](http://xkcd.com/221/)
[![Tour of Accounting](http://assets.amuniversal.com/321a39e06d6401301d80001dd8b71c47)](http://dilbert.com/strip/2001-10-25)
