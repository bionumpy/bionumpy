Position Weight Matrix
-------------------------

This tutorial shows how to use rollable functions. Reading a motif-pwm from file, a PositionWeightMatrix function is created using the appropriate alphabet and counts. Since `PositionWeightMatrix` is a `RollableFunction` subclass it has a rolling_window method that applies the pwm to all valid windows in the sequence set.

.. literalinclude:: /../scripts/pwm_example.py
