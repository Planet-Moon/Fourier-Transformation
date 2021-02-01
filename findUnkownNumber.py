findNumber = 1e8
currentNumber = 1
lesser = 0
larger = 0
numberFound = False
numbersTested = []
actionsTested = []
between = []

while not numberFound:
    numbersTested.append(currentNumber)
    if abs(currentNumber - findNumber) < 0.001:
        numberFound = True
        actionsTested.append("found")
    elif currentNumber < findNumber:
        lesser = currentNumber
        if larger < lesser:
            currentNumber *= 2
            actionsTested.append("double")
        else:
            currentNumber = (lesser + larger)/2
            actionsTested.append("half more")
    elif currentNumber > findNumber:
        larger = currentNumber
        currentNumber = (lesser + larger)/2
        actionsTested.append("half less")
    between.append([lesser, larger])
    pass
    
print("Number found in {} steps".format(len(numbersTested)))
for i in range(len(actionsTested)):
    print("step {}: number: {}\t{}\tbetween: {} and {}" 
        .format(
            i+1,
            numbersTested[i],
            actionsTested[i],
            between[i][0],
            between[i][1]
        )
    )
    