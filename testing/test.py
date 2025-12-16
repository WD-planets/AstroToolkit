from ATK.utilities.defaults import RETURNS

temp = RETURNS.EXCEPTION


if temp is RETURNS.EXCEPTION:
    print("temp IS RETURNS.EXCEPTION")

if temp in (RETURNS.EXCEPTION, RETURNS.NULL):
    print("temp IS IN (RETURNS.EXCEPTION, RETURNS.NULL)")

if temp == RETURNS.EXCEPTION:
    print("temp == RETURNS.EXCEPTION")
