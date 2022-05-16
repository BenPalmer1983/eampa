#!/bin/python3

import numpy
import time
#from f_toolbox import atom

from g import g
from ds import ds



class std:



  @staticmethod
  def read(file_path):

    #   Read file
    ###############################################################

    block = ''
    fh = open(file_path, 'r')
    for line in fh:
      block = block + line
    fh.close()

    no_comments = ''
    in_comment = 0
    n = 0
    while(n<len(block)):
      if(in_comment == 0 and n < len(block) - 1 and block[n:n+2] == "/*"):
        in_comment = 1
      elif(in_comment == 0 and n < len(block) - 1 and block[n:n+2] == "//"):
        in_comment = 2
      elif(in_comment == 1 and n < len(block) - 1 and block[n:n+2] == "*/"):
        in_comment = 0
        n = n + 1
      elif(in_comment == 2 and n < len(block) and block[n] == "\n"):
        in_comment = 0
        no_comments = no_comments + block[n]
      elif(in_comment == 0):
        no_comments = no_comments + block[n]
      n = n + 1

    d = []
    lines = no_comments.split('\n')
    for line in lines:
      line = std.one_space((line.replace("\t", " ")).strip())
      if(line != ""):
        d.append(line)  
    return d

  @staticmethod
  def one_space(line, sep=" "):
    out = ''   
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0 and not (char == " " and last_char == " ")):
        out = out + char
      last_char = char
    return out 














  """
  @staticmethod
  def one_space(inp):
    out = ''
    last = None
    for c in inp:
      if(not (last == " " and c == " ")):
        out = out + c
      last = c
    return out






  """