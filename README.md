# screenr
An R package to enable screening in/out subjects who are likely to to
test positive/negative, respectively.

Consider the situation where a definitive test for some condition
is expensive, and the condition is rare.  In that case, universal testing
would not be efficient in terms of the yield of postive results per test
performed.  Now suppose that responses to a set of simple diagnostic
questions or observations may be predictive of the definitive test result.
Package screenr enables the estimation of thresholds for making decisions
about when to perform the definitive test on newly observed subjects based
on Receiver Operating Characteristics (ROC) estimated from an initial sample.
The choice of a particular screening threshold is left to the user, and
should be based on careful consideration of application-specific tradeoffs
between sensitivity (true positive fraction) and specificity (true negative
fraction).

See also: Tutorial_Example.html
