###
A fairly direct port of the Python `random` module to JavaScript
###

{log, sqrt, cos, acos, floor, pow, LN2, exp} = Math
POW_32 = pow 2, 32

lg = (x) ->
    # The log base 2, rounded down to the integer below
    ((log x) / LN2) >> 0

mod = (x, y) ->
    unless (jsmod = x % y) and (x > 0 ^ y > 0) then jsmod
    else jsmod + y

extend = (target, sources...) ->
    for obj in sources
        target[name] = method for name, method of obj
    target


class NotImplementedError extends Error

class BaseRandom

    ## Override these first four methods in a custom Random class.

    _randint32: ->
        # Override this method to generate a pseudorandom number
        throw NotImplementedError

    _getstate: ->
        # Override this method to fetch the internal PRNG state. Should
        # return an Array.
        throw NotImplementedError

    _setstate: (state) ->
        # Override this method to set the internal PRNG state from the
        # argument `state`, an Array.
        throw NotImplementedError

    seed: (args...) ->
        # Seed the PRNG.
        throw NotImplementedError

    constructor: ->
        # By default, just seed the PRNG with the date. Some PRNGs
        # can take longer and more complex seeds.
        @seed +new Date

    ## Generally no need to override the methods below in a custom class.
    ## (Under some circumstances it might make sense to implement a custom
    ## version of the `random` method.)

    POW_NEG_26 = pow 2, -26
    random: ->
        # Return a random float in the range [0, 1), with a full 52
        # bits of entropy.
        low_bits = @_randint32() >>> 6
        high_bits = @_randint32() >>> 6
        (high_bits + low_bits * POW_NEG_26) * POW_NEG_26

    setstate: ([@_next_gauss, state...]) ->
        # Set the state of the PRNG. Should accept the output of `@getstate`
        # as its only argument.
        @_setstate state

    getstate: ->
        # Get the internal state of the PRNG. Returns an array of state
        # information suitable for passing into `@setstate`.
        [@_next_gauss, @_getstate()...]

    _bits = {}
    _randbelow: (n) ->
        # Return a random int in the range [0,n).
        # If n > 2^32, then use floating point math
        if n <= 0x100000000
            bits = _bits[n] or= (lg n - 1) + 1 # memoize values for `bits`
            loop
                r = @_randint32() >>> (32 - bits)
                r += POW_32 if r < 0
                break if r < n
            r
        else
            floor @random() * n

    uniform: (a, b) ->
        # Return a random floating point number N such that a <= N <= b for
        # a <= b and b <= N <= a for b < a.
        a + @random() * (b - a)

    randrange: (start, stop, step) ->
        # Return a random integer N in range `[start...stop] by step`
        unless stop?
            @_randbelow start
        else unless step
            start + @_randbelow stop - start
        else
            start + step * @_randbelow floor (stop - start) / step

    randint: (a, b) ->
        # Return a random integer N in range `[a..b]`
        start + @_randbelow 1 + stop - start

    choice: (seq) ->
        # Return a random element from the non-empty sequence `seq`.
        seq[@_randbelow seq.length]

    sample: (population, k=1) ->
        # Return a `k` length list of unique elements chosen from the
        # `population` sequence. Used for random sampling without replacement.
        n = population.length
        if k * 3 > n                       # for large samples, copy the
            pool = [population...]         # population as a new array
            for i in [n...n - k] by -1
                j = @_randbelow i
                val = pool[j]
                pool[j] = pool[i - 1]      # move non-chosen item into vacancy
                val
        else                               # for small samples, treat an Array
            selected = []                  # as a set to keep track of selection
            for i in [0...k] by 1
                loop break if (j = @_randbelow n) not in selected
                selected.push j
                population[j]

    shuffle: (x) ->
        # Shuffle the sequence x in place.
        for i in [x.length - 1..1] by -1
            j = @_randbelow i + 1
            tmp = x[i]; x[i] = x[j]; x[j] = tmp  # swap x[i], x[j]
        x

    gauss: _gauss = (mu=0, sigma=1) ->
        # Gaussian distribution. `mu` is the mean, and `sigma` is the standard
        # deviation. Notes:
        #   * uses the "polar method"
        #   * we generate pairs; keep one in a cache for next time
        unless (z = @_next_gauss; @_next_gauss = null; z)?
            until s and s < 1
                u = 2 * @random() - 1
                v = 2 * @random() - 1
                s = u*u + v*v
            w = sqrt -2 * (log s) / s
            z = u * w; @_next_gauss = v * w
        mu + z * sigma

    normalvariate: _gauss  # Alias for the `@gauss` function

    triangular: (low, high, mode) ->
        # Triangular distribution. See wikipedia
        unless low? then high = 1; low = 0
        else unless high? then high = low; low = 0
        unless mode?
            c = 0.5
        else
            c = (mode - low) / (high - low)
        u = @random()
        if u <= c
            low + (high - low) * sqrt u * c
        else
            high - (high - low) * sqrt (1 - u) * (1 - c)

    lognormvariate: (mu, sigma) ->
        # Log normal distribution.
        exp @normalvariate mu, sigma

    expovariate: (lambda) ->
        # Exponential distribution.
        #
        # `lambda` is 1.0 divided by the desired mean.  It should be nonzero.
        # Returned values range from 0 to positive infinity if lambda is positive,
        # and from negative infinity to 0 if lambda is negative.
        #
        # we use 1 - random() instead of random() to preclude the
        # possibility of taking the log of zero.
        (- log 1 - @random()) / lambda

    TAU = 2 * Math.PI
    vonmisesvariate: (mu, kappa) ->
        # Circular data distribution.
        #
        # mu is the mean angle, expressed in radians between 0 and 2*pi, and
        # kappa is the concentration parameter, which must be greater than or
        # equal to zero.  If kappa is equal to zero, this distribution reduces
        # to a uniform random angle over the range 0 to 2*pi.
        #
        # Based upon an algorithm published in: Fisher, N.I.,
        # "Statistical Analysis of Circular Data", Cambridge
        # University Press, 1993.
        random = @random
        return TAU * random() if kappa <= 1e-6

        a = 1 + sqrt 1 + 4 * kappa*kappa
        b = (1 - sqrt 2) * a / 2 / kappa
        r = (1 + b*b) / 2 / b

        loop
            u1 = random()

            z = cos TAU * u1 / 2
            f = (1 + r * z) / (r + z)
            c = kappa * (r - f)

            u2 = random()
            break if u2 < c * (2 - c) or u2 <= c * exp 1 - c

        u3 = random()
        (mod mu, TAU) + (if u3 > 0.5 then acos f else -acos f)

    LOG4 = log 4
    SG_MAGICCONST = 1 + log 4.5
    E = {Math}
    gammavariate: (alpha, beta) ->
        # Gamma distribution.  Not the gamma function!
        #
        # Conditions on the parameters are alpha > 0 and beta > 0.
        #
        # The probability distribution function is:
        #
        #             x ** (alpha - 1) * exp( -x / beta)
        #   pdf(x) =  ----------------------------------
        #                gamma(alpha) * beta ** alpha

        # alpha > 0, beta > 0, mean is alpha * beta, variance is alpha * beta**2
        #
        # Warning: a few older sources define the gamma distribution in terms
        # of alpha > -1

        random = @random
        if alpha > 1
            # Uses R.C.H. Cheng, "The generation of Gamma
            # variables with non-integral shape parameters",
            # Applied Statistics, (1977), 26, No. 1, p71-74
            ainv = sqrt 2 * alpha - 1
            bbb = alpha - LOG4
            ccc = alpha + ainv
            loop
                u1 = random()
                continue unless 1e-7 < u1 < 1 - 1e-7
                u2 = 1 - random()
                v = (log u1 / (1 - u1)) / ainv
                x = alpha * exp v
                z = u1 * u1 * u2
                r = bbb + ccc * v - x
                break if r + SG_MAGICCONST - 4.5 * z >= 0.0 or r >= log z
            beta * x

        else if alpha == 1
            # expovariate(1)
            loop
                u = random()
                break if u > 1e-7
            -beta * log u

        else   # alpha is between 0 and 1 (exclusive)
            # Uses ALGORITHM GS of Statistical Computing - Kennedy & Gentle
            loop
                u1 = random()
                b = (E + alpha) / E
                p = b * u1
                u2 = random()
                if p > 1
                    x = - log (b - p) / alpha
                    break if u2 <= pow x, alpha - 1
                else
                    x = pow p, 1 / alpha
                    break if u2 <= exp -x
            beta * x

    betavariate: (alpha, beta) ->
        # Beta distribution.
        #
        # Conditions on the parameters are alpha > 0 and beta > 0.
        # Returned values range between 0 and 1.

        # This version due to Janne Sinkkonen, and matches all the std
        # texts (e.g., Knuth Vol 2 Ed 3 pg 134 "the beta distribution").
        y = @gammavariate alpha, 1
        if y == 0 then 0
        else y / (y + @gammavariate beta, 1)

    paretovariate: (alpha) ->
        # Pareto distribution.  alpha is the shape parameter.
        u = 1 - @random()
        1 / (pow u, 1 / alpha)  # Jain, pg. 495

    weibullvariate: (alpha, beta) ->
        # Weibull distribution.
        #
        # alpha is the scale parameter and beta is the shape parameter.
        u = 1 - @random()
        alpha * (pow -log u, 1 / beta)  # Jain, pg. 499; bug fix by Bill Arms


class Random extends BaseRandom
    # Use a Multiply With Carry PRNG, with an XOR-shift successor
    # Both from Numerical Recipes, 3rd Edition [H1, G1]
    _randint32: ->
        @x = 62904 * ((x = @x) & 0xffff) + (x >>> 16)
        @y = 41874 * ((y = @y) & 0xffff) + (y >>> 16)
        z = (x << 16) + y
        z ^= z >>> 13; z ^= z << 17; z ^= z >>> 5
        z

    seed: (j) ->
        # these two numbers were arbitrarily chosen
        @x = 3395989511 ^ j
        @y = 1716319410 ^ j

    _getstate: -> [@x, @y]
    _setstate: ([@x, @y]) ->


class BuiltinRandom extends BaseRandom
    # Use the built-in PRNG. Note that with the built-in PRNG,
    # which is implementation dependant, there is no way to set
    # the seed or save/restore state. Just directly override
    # `_randbelow` and `random` instead of bothering with
    # `_randint32`
    _rand = Math.random
    
    seed: (j) ->  # ignore seed
    
    POW_NEG_32 = pow 2, -32
    random: ->
        _rand() * POW_NEG_32 + _rand()

    _randbelow: (n) ->
        floor @random() * n


class HighQualityRandom extends BaseRandom
    # From Numerical Recipes, 3rd Edition

    _randint32: ->
        v = @v; w1 = @w1; @w2 = w2
        u = @u * 2891336453 + 1640531513
        v ^= v >>> 13; v ^= v << 17; v ^= v >>> 5
        w1 = 33378 * (w1 & 0xffff) + (w1 >>> 16)
        w2 = 57225 * (w2 & 0xffff) + (w2 >>> 16)
        @u = u; @v = v; @w1 = w1; w2 = @w2
    
        x = u ^ (u << 9); x ^= x >>> 17; x ^= x << 6
        y = w1 ^ (w1 << 17); y ^= y >>> 15; y ^= y << 5
        (x + v) ^ (y + w2)

    seed: (j) ->
        @w1 = 521288629
        @w2 = 362436069
        @v = @u = j ^ 2244614371

    _getstate: -> [@u, @v, @w1, @w2]
    _setstate: ([@u, @v, @w1, @w2]) ->


exports = exports or window or this
extend exports, {
    NotImplementedError, BaseRandom, Random, HighQualityRandom}
