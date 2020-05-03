# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 13:03:17 2018

@author: aesseb
"""

#I'm a comment...

####### Data Types #######
## Not all things are made equal, operations can be done to a string that
## cannot be done to an integer

#int (integer)
5
type(5)
#float
5.0
#Boolean
True
#String
"Stuff"
"Characters $*@"
type("Hello")
#List
[0,1,2]
#Dictionary
{"Bioinformatics":"awesome", "Politics":"meh"}

####### Variable Assignment #######
##Creating a variable is giving a data type (aka an object) a name ##
x = 5
y = "Hello"
i_like_bio = True
a = x
my_list = [1,2,3,4,5]

####### Print Statements #######
## To display results, help track how variables change and to debug code
print(x)
print(y)
print(a)
print(y, a, i_like_bio)

#### COMPLETE EXERCISE 1 IN THE GUIDE TO REINFORCE THESE CONCEPTS ####

####### Flow Control #######

####### Conditional statements #######
# if condition or variable is True:
    # execute code written here
# else:
    # condition is False, execute code written here
    
#Example conditions:
if 10 == 11:
    print (True)
else:
    print (False)
    
if 'a' == 'ab':
    print (True)
else:
    print ("These strings are not equal")
    
if 1 in [1,2,3]:
    print ("Item in list")
#You do not always have to include an else statement

if True # This line has an error because the colon is missing
    print (50)

if apple: # This line will cause an error because the variable apple does not exist
    print ("Pear")
    
apple = True
if apple:
    print ("Pear") # Now it works because apple has been set to True
    
test = 80
if test < 10 :
    print ("Less than 10")
elif test < 50 : #When you have more than two conditions you require the elif statement
    print ("Less than 50")
else:
    print ("Greater than or equal to 50")
    
#### COMPLETE EXERCISE 2 IN THE GUIDE TO REINFORCE THESE CONCEPTS ####

####### Lists #######
dnalist = ["A", "T", "G", "C"] #Square brackets define a list
print (len(dnalist))
#Below, we use square brackets to extract what's stored at the specified position in the list
#Python ALWAYS counts from 0 (not 1)
# This is known as indexing
print (dnalist[0], dnalist[3]) 
print (dnalist[1:3])

# Strings are like a list of characters so, you can do similar things with strings
dnastr = "ATGC"
print (len(dnastr))
print (dnastr[0], dnastr[3]) 
print (dnastr[1:3])

####### Loops #######

i = 0
while i < 8:
    print ("Grade = ", i)
    #If you do not change the value of i, the loop will continue forever, or until your computer crashes!
    i += 1 #Increases value of i by 1

for i in range(8):
    print ("Grade = ", i)
    
""" The two loops above have the exact same outcome but look different in terms of how they're written.
     You can use loops to ITERATE over lists, dictionaries, strings and other data structures 
     which are iterable """
    
####### Exercise - try to figure out this code #######

import random # A Python module that can generate random numbers
# When you import a module, you can then access all the code/functions/classes that
# are contained within the module

for i in range(10):
    grade = random.randint(1,7)
    if grade == 7:
        print ("Very happy")
    elif grade >= 4:
        print ("I passed") 
    else:
        print ("Not so happy")

# Try adding a print statement to the above for loop to determine how the variable
# grade changes

######## How to build lists using loops ########

mylist = [] #An empty list
for i in range(10):
    random_number = random.randint(0, 3)
    mylist.append(dnalist[random_number]) # Append adds item in brackets to the end of a mylist
print (mylist)

# If you run the above code multiple times, the answer will change. Can you figure out why?

####### Dictionaries #######
# Have a key:value relationship with unique keys

transfer = {"lan": "fast", "wi-fi":"fine", "4g":"OK", "key":"value"} #Dicsionaries are defined using curly brackets
approach = "wi-fi"
print (approach, "transfer is", transfer[approach]) #Use square brackets and the key to access a value
print("But lan is", transfer["lan"])

#Empty dictionary
emptyDict = {}
#Add key:value pair
emptyDict["key"] = "value"
emptyDict["A"] = 1
emptyDict["B"] = 2

# How to loop through items in dictionary
for key, value in transfer.items():
    print(key, value)
    
####### Functions #######
# A function is a piece of code that can be called  and used in many different places
# Instead of writing the same piece of code many times, you write it once and make it a function
# It’s like a variable name for a whole piece of  code rather than one value

def multiply(a, b):
    return a*b # terminates the function and passes this value back to where the function was called

def multiply_2(a, b):
    print(a*b)
    # Without a return statement, nothing will be passed back to where this function was called
    # Therefore, below, b will not have a value and the result will be printed out but not stored anywhere

def multiply_3(a, b, c):
    print(a*b*c)
    return a*b*c # This function will both print the result and return the result to where it was called

a = multiply(4, 6)
print (a)

b = multiply_2(3, 9)
print (b)

c = multiply_3(10, 50, 29)
print (c)

#### TRY EXERCISE 4 IN THE GUIDE ####

####### Classes #######
#Should be read about in the guide with guide.py providing numerous examples of classes




      

    


