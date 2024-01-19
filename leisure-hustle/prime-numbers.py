# Returns a list of prime numbers #
#############################################################
from IPython.display import clear_output
'''
definitionwise, a prime number is only divisible by 
1 and the number itself. that means, no other number 
can divide it with a remainder of zero. so, as every
number between 1 and selected number tries to divide
the selected number, the remainder would never equal
zero.

that means, the number of divisors who yield remainder 
unequal to zero would be the total number of divisors.

first, take all the numbers between 1 and selected 
number inside a list. second, run a loop and make sure when
each number yields remainder yields non-zero remainder, 
they are also appended in another list. 

compare the lists and if the lengths of the lists turn
out to be equal, then the selected number must be a prime
number. 
'''
# initialize a list to store all prime numbers
primelist = []

# get the input from the user
num = int(input('What range do you want, pal?'))
clear_output(wait=True)

# run the loop
for i in range(1, num+1):
  nondivisors = []
  divisors = list(range(2, i))
  for d in divisors:
    if i%d != 0:
      nondivisors.append(d)
      
  if len(divisors) == len(nondivisors):
    primelist.append(i)

primelist





# Returns if a number is prime or not #
#############################################################
from IPython.display import clear_output

# get the number input from the user
num = int(input('Hey pal, tell me which number to check for prime?'))
clear_output(wait=True)

# run the condition loop
for i in range(1, num+1):
  nondivisors = []
  divisors = list(range(2, i))
  for d in divisors:
    if i%d != 0:
      nondivisors.append(d)

# check if the condition was met
if len(divisors) == len(nondivisors):
  print(f'{num}: Prime')
else:
  print(f'{num}: Non-prime')
