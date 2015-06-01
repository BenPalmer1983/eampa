# !/usr/bin/python
import os, sys, time

def trim(inputStr):
  output = ""
  i = len(inputStr)
  lastSpace = 0
  while i>0:
    charSelected = subString(inputStr,i,1)
    if lastSpace==0:
      if charSelected!=" ":
        lastSpace = 1
        output = charSelected + output
    else:    
      output = charSelected + output
    i = i - 1
  return output
def removeCR(inputStr):
# Removes spaces
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) not in (10,13):
      output = output + testChar
  return output     
def removeTabs(inputStr):
# Removes spaces
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) != 9:
      output = output + testChar
  return output   
def trimLeading(inputStr):
  output = ""
  i = 0
  iEnd = len(inputStr)
  firstSpace = 1
  while i<iEnd:
    i = i + 1
    charSelected = subString(inputStr,i,1)
    if firstSpace==1:
      if ord(charSelected)!=(32):
        firstSpace = 0
        output = output + charSelected
    else:    
      output = output + charSelected
  return output      
def spacesString(length):
  i = 0
  output = ""
  while i<length:  
    i = i + 1
    output = output + " "
  return output      
def subString(inputStr, start, length):
  a = start-1
  output = inputStr[a:a+length]
  return output
  
#input arguments
arguments=[None]*10
counter = 0
for arg in sys.argv:
  counter = counter + 1
  arguments[counter] = arg
  
# 1 = remove tabs etc, format the layout of the code, add in compiled time stamp  
# 2 = check compiled OK
# 3 = Remove unused variables


#========================================================================
# Format layout, remove unnecessary tabs, spaces
#========================================================================

if arguments[2]=="1":
 print (time.strftime("%H:%M:%S  %d/%m/%Y"))
 # Loop through files, remove tabs, trailing spaces etc
 directory = "./src"
 for file in os.listdir(directory):
  spaces = 0
  inFunction = 0 # Inside function or subroutine
  if file.endswith(".f90"):
    srcFile = open(directory+"/"+file, 'r')
    outputFile = open(directory+"/"+file+".temp", 'w')
    lastLine = ""
    for line in srcFile: 
      newLine = line      
# Remove cr, tabs and trailing spaces    
      newLine = removeCR(newLine)  # Remove carriage returns
      newLine = removeTabs(newLine) # Remove tabs
      newLine = trim(newLine)      # Remove trailing spaces   
      inputLine = newLine
# change common statements   
      newLine = trimLeading(newLine)
      newLineU = newLine.upper()
      lineLen = len(newLine)
      lineMod = 0
# Save program name      
      if newLineU=="PROGRAM" or subString(newLineU,1,8)=="PROGRAM ": 
        lineArray = newLine.split(" ")  
        print lineArray[0],lineArray[1]
        programName = lineArray[1]       
        newLine = "PROGRAM"+subString(newLine,8,lineLen-7) 
        lineMod = 1
        lineLen = len(newLine)
# Save program name (END)     
      if(newLineU=="END"):
        newLine = "End Program "+programName
        lineMod = 1
        lineLen = len(newLine)
# Update compile line in globals.f90, if exists
      if(file=="globals.f90"):
        if(subString(newLine,1,11)=="compileLine"):
          newLine = "compileLine = "+chr(34)+time.strftime("%H:%M:%S  %d/%m/%Y")+chr(34)
          lineMod = 1
          lineLen = len(newLine)
# If      
      if subString(newLineU,1,3)=="IF(":
        newLine = "If("+subString(newLine,4,lineLen-3)
        lineMod = 1
        lineLen = len(newLine)
      elif subString(newLineU,1,4)=="IF (":
        newLine = "If("+subString(newLine,5,lineLen-4)
        lineMod = 1
        lineLen = len(newLine)
# Else      
      if subString(newLineU,1,4)=="ELSE":
        newLine = "Else"+subString(newLine,5,lineLen-4)
        lineMod = 1
        lineLen = len(newLine)
# Else If      
      if subString(newLineU,1,9)=="ELSE IF (":
        newLine = "ElseIf("+subString(newLine,10,lineLen-9)
        lineMod = 1
        lineLen = len(newLine)
      elif subString(newLineU,1,8)=="ELSE IF(":
        newLine = "ElseIf("+subString(newLine,9,lineLen-8)
        lineMod = 1
        lineLen = len(newLine)
      elif subString(newLineU,1,8)=="ELSEIF (":
        newLine = "ElseIf("+subString(newLine,9,lineLen-8)
        lineMod = 1
        lineLen = len(newLine)
      elif subString(newLineU,1,7)=="ELSEIF(":
        newLine = "ElseIf("+subString(newLine,8,lineLen-7)
        lineMod = 1
        lineLen = len(newLine)
# End If      
      if subString(newLineU,1,6)=="END IF":
        newLine = "End If"+subString(newLine,7,lineLen-6)
        lineMod = 1
        lineLen = len(newLine)
      elif subString(newLineU,1,5)=="ENDIF":
        newLine = "End If"+subString(newLine,6,lineLen-5)
        lineMod = 1
        lineLen = len(newLine)
# Do     
      if subString(newLineU,1,3)=="DO ":
        newLine = "Do "+subString(newLine,4,lineLen-3)
        lineMod = 1
        lineLen = len(newLine)
# End Do     
      if subString(newLineU,1,6)=="END DO":
        newLine = "End Do"+subString(newLine,7,lineLen-6)
        lineMod = 1
        lineLen = len(newLine)
      elif subString(newLineU,1,5)=="ENDDO":
        newLine = "End Do"+subString(newLine,6,lineLen-5)
        lineMod = 1
        lineLen = len(newLine)
# Then
      if subString(newLineU,lineLen-3,4)=="THEN": 
        newLine = subString(newLine,1,lineLen-4)+"Then"
        lineMod = 1  
        lineLen = len(newLine)   
      if subString(newLineU,lineLen-4,5)==" THEN": 
        newLine = subString(newLine,1,lineLen-5)+"Then"
        lineMod = 1 
        lineLen = len(newLine)        
# Module     
      if subString(newLineU,1,6)=="MODULE":
        newLine = "Module"+subString(newLine,7,lineLen-6)
        lineMod = 1   
        lineLen = len(newLine)       
# Subroutine     
      if subString(newLineU,1,10)=="SUBROUTINE":
        newLine = "Subroutine"+subString(newLine,11,lineLen-10)
        lineMod = 1     
        inFunction = 1  
        lineLen = len(newLine)
# End Subroutine     
      if subString(newLineU,1,14)=="END SUBROUTINE":
        newLine = "End Subroutine"+subString(newLine,15,lineLen-14)
        lineMod = 1     
        inFunction = 0 
        lineLen = len(newLine) 
      elif subString(newLineU,1,13)=="ENDSUBROUTINE":
        newLine = "End Subroutine"+subString(newLine,14,lineLen-13)
        lineMod = 1  
        inFunction = 0   
        lineLen = len(newLine)
# Function     
      if subString(newLineU,1,8)=="FUNCTION":
        newLine = "Function"+subString(newLine,9,lineLen-8)
        lineMod = 1   
        inFunction = 1   
        lineLen = len(newLine) 
# End Function     
      if subString(newLineU,1,12)=="END FUNCTION":
        newLine = "End Function"+subString(newLine,13,lineLen-12)
        lineMod = 1       
        inFunction = 0  
        lineLen = len(newLine)
# Recursive Function     
      if subString(newLineU,1,18)=="RECURSIVE FUNCTION":
        newLine = "Recursive Function"+subString(newLine,19,lineLen-18)
        lineMod = 1     
        inFunction = 1  
        lineLen = len(newLine) 
# Update uppercase
      newLineU = newLine.upper()
# Calculate new leading spaces (reduce spaces)
      if subString(newLineU,1,4)=="END ": 
        spaces = spaces - 2 
      if subString(newLineU,1,4)=="ELSE":  
        spaces = spaces - 2        
# Replace leading spaces        
      if(subString(newLineU,1,1)!="!"):
        newLine = spacesString(spaces)+newLine
      else:
        if(subString(newLineU,1,2)!="! "):
          newLine = "! "+subString(newLine,2,lineLen-1)
      newLine = trim(newLine)   
# Calculate new leading spaces (increase spaces)     
      if subString(newLineU,1,7)=="PROGRAM": 
        spaces = spaces + 2 
      if subString(newLineU,1,6)=="MODULE": 
        spaces = spaces + 2 
      if subString(newLineU,1,10)=="SUBROUTINE":
        spaces = spaces + 2
      if subString(newLineU,1,3)=="IF(":  
        spaces = spaces + 2
      if subString(newLineU,1,4)=="ELSE":  
        spaces = spaces + 2
      if subString(newLineU,1,3)=="DO ":  
        spaces = spaces + 2
      if subString(newLineU,1,9)=="FUNCTION ":  
        spaces = spaces + 2
      if subString(newLineU,1,19)=="RECURSIVE FUNCTION ":  
        spaces = spaces + 2
# write line
      if inFunction==1:
        if newLine!="":
          outputFile.write(newLine+"\n")
      if inFunction==0:
        if lastLine!="":
          outputFile.write(newLine+"\n")
        else:  
          if newLine!="":
            outputFile.write(newLine+"\n")
#store "last loop" line
      lastLine = newLine    
# End loop      
    os.remove(directory+"/"+file)
    os.rename(directory+"/"+file+".temp",directory+"/"+file)     
    
#========================================================================
# Check compiled OK, or with fatal error
#========================================================================    

if arguments[2]=="2":    
  errorFile = open('./logs/make.log', 'r')
  compileStatus = 1
  for line in errorFile:
    newLine = trim(line.upper())
    if 'FATAL ERROR' in newLine:
      compileStatus = 0
  
  if compileStatus==1:
    print "Compiled OK"
  else:  
    print "Compiled with Fatal Error"
    
 
#========================================================================
# Remove unused variables
#========================================================================

if arguments[2]=="3":
    #print "Reading File"
    errorFile = open('./logs/build.log', 'r')
    fileRows=[None]*10000
    counter = 0
    compileStatus = 1
    for line in errorFile:
      errLine = trim(line.upper())
      if 'FATAL ERROR' in errLine:
        compileStatus = 0
      rowStripped = line.strip()
      if rowStripped:
        counter = counter + 1
        fileRows[counter] = rowStripped
    fileRowCount = counter
    if compileStatus==1:
        print "Compiled OK"
    else:  
        print "Compiled with Fatal Error"
    fileListArr=[None]*10000
    fileNameArr=[None]*10000
    varRowArr=[None]*10000
    varNameArr=[None]*10000
    counter = 0
    counterFN = 0
    for j in range(1,fileRowCount+1):
      if 'Unused variable' in fileRows[j]:
        counter = counter + 1
        line = fileRows[j]
        strRemainder,variableName,strRemainder = line.split("'",2)
        line = fileRows[j-3]    
        fileName,varDetails,strRemainder = line.split(":",2)
        varLine, varPos = varDetails.split(".",2)    
    # File name
        if counterFN==0:
          counterFN = 1
          fileListArr[counterFN] = "./"+fileName
          print fileListArr[counterFN]
        else:
          if fileListArr[counterFN]!="./"+fileName:
            counterFN = counterFN + 1
            fileListArr[counterFN] = "./"+fileName
    # Variable array        
        fileNameArr[counter] = "./"+fileName
        varRowArr[counter] = varLine
        varNameArr[counter] = variableName
    #    print "Unused ",variableName
    varCount = counter
    fileCount = counterFN
    for j in range(1,fileCount+1):                    # Loop through files
      srcFile = open(fileListArr[j], 'r')
      outputFile = open(fileListArr[j]+".temp", 'w')
      lineCounter = 0 
      for line in srcFile:                          # Loop through lines in the file  
        lineCounter = lineCounter + 1
        inputLine = trimAll(line)
    # Check line for amendment
        amendment = 0
        for i in range(1,varCount+1):                 # Loop through vars to ammend/remove
          if fileNameArr[i]==fileListArr[j]:
            if int(varRowArr[i])==int(lineCounter): 
    # split line, list of variables
              varToRemove = varNameArr[i]
              varType, varList = inputLine.split("::",2)
              varList = removeBlanks(varList)
              varArr = varList.split(",")
    # loop through variables and check for var to remove
              varListNew = ""
              for k in range(0,len(varArr)):
                varNameUpper = removeBlanks(varToRemove.upper())
                varNameTestUpper = removeBlanks(varArr[k].upper())
                if varNameUpper!=varNameTestUpper:
                  varListNew = varListNew + trimAll(varArr[k]) + ","
              varListNew = trimAll(varListNew)   
              if varListNew:
                newLine = trimAll(varType)+" :: "+varListNew   
              else:  
                newLine = ""  
    # remove trailing comma  
              newLine = trimAll(newLine)            
              newLine = removeTrailingComma(newLine)    
              amendment = 1
    # store as input line, for other same line amendments
              if newLine: 
                inputLine = newLine
        if amendment==0:      
          #print trimAll(line)
          outputFile.write(inputLine+"\n")
        else:      
          outputFile.write(newLine+"\n")
    #Replace src files with new files    
    for j in range(1,fileCount+1): 
      os.remove(fileListArr[j])
      os.rename(fileListArr[j]+".temp",fileListArr[j])

   
    