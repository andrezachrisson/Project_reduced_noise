#!/usr/bin/env python3
def noise_checker(column):
    if indel_counter(column) and uniq(column) and amino_acid_twice(column):
        return True
    else:
        return False

def indel_counter(column):
    indels = column.count('-')
    column_length = len(column)
    value = indels/column_length
    if value < 0.5:
        return True
    else:
        return False

def uniq(column):
    dict = {}
    count = 0
    for amino in column:
        if amino != "-" and amino not in dict:
            count +=1
            dict[amino] = count
    try:
        list(filter(lambda a: a != "-", column))
    except:
        pass

    if count/len(column) > 0.5:
        return False
    else:
        return True

def amino_acid_twice(column):
    dict = {}
    for i in column:
        if i != '-' and i not in dict:
            value = column.count(i)
            dict[i] = value
            if value > 2:
                return True
                break
    return False
