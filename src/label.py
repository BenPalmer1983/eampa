#!/bin/python3

import numpy

from g import g
from ds import ds


class label:

  @staticmethod
  def set(inp):
    inp = inp.strip().upper()
    for n in range(len(g.labels)):
      if(g.labels[n] == inp):
        return n
    g.labels.append(inp)
    return (len(g.labels) - 1)

  @staticmethod
  def get(inp):
    return g.labels[inp]
