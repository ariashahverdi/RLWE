#!/usr/bin/python

# Parses RLWE challenge/instance files:
# takes the path to a .challenge file first, and an optional second path to a
# corresponding .instance file. Prints the contents of the .challenge file
# and the .instance file, if a path is provided.
# See the code for examples of how to access message contents.
#
# You will need the "protobuf" python package to run this code, which you can
# install with "easy_install protobuf"

# You can find more detailed information about parsing messages in python here:
# https://developers.google.com/protocol-buffers/docs/reference/python-generated
import sys
import Challenges_pb2
import math


def parse_challenge(challenge_path):
    if challenge_path is None:
        return
    challenge = None
    with open(challenge_path) as f:
        challenge = Challenges_pb2.Challenge()
        challenge.ParseFromString(f.read())

    return challenge


def parse_instance(challenge, instance_path):
    if challenge is None or instance_path is None:
        return

    instance = None

    with open(instance_path) as f:
        paramType = challenge.WhichOneof("params")
        if paramType == "cparams":
            instance = Challenges_pb2.InstanceCont()
        elif paramType == "dparams":
            instance = Challenges_pb2.InstanceDisc()
        elif paramType == "rparams":
            instance = Challenges_pb2.InstanceRLWR()
        else:
            return
        instance.ParseFromString(f.read())
    return instance


def parse_secret(secret_path):
    if secret_path is None:
        return
    secret = None
    with open(secret_path) as f:
        secret = Challenges_pb2.Secret()
        secret.ParseFromString(f.read())
    return secret


def reverseBits(num, bitSize):
    # convert number into binary representation
    # output will be like bin(10) = '0b10101'
    binary = bin(num)

    # skip first two characters of binary
    # representation string and reverse
    # remaining string and then append zeros
    # after it. binary[-1:1:-1]  --> start
    # from last character and reverse it until
    # second last character from left
    reverse = binary[-1:1:-1]
    reverse = reverse + (bitSize - len(reverse)) * '0'

    # converts reversed binary string into integer
    # print int(reverse,2)
    return int(reverse, 2)


def adjust_index(x):
    res = [0 for _ in range(len(x))]
    bitSize = int(math.log(len(x), 2))
    for i in range(len(x)):
        res[reverseBits(i, bitSize)] = x[i]
    return res


def get_challenge(challenge_path, instance_path, secret_path):

    challenge = parse_challenge(challenge_path)

    if challenge is None:
        print "Could not parse the challenge parameters."
        sys.exit(-1)
    print challenge # print all contents

    # examples of how to access message contents
    # number of instances in the challenge
    #print challenge.numInstances
    # shows which type of parameter the message contains
    #print challenge.WhichOneof("params")
    # cyclotomic index of a challenge
    #print {'cparams': challenge.cparams.m,
    #       # number of samples in each instance
    #       'dparams': challenge.dparams.numSamples,
    #       # input modulus of the RLWR challenge
    #       'rparams': challenge.rparams.q}[challenge.WhichOneof("params")]

    secret = parse_secret(secret_path)

    if secret is None:
        print "Could not parse secret."
        sys.exit(-1)
    else:
        s = adjust_index(secret.s.xs[0:])
        #print secret  # print all contents
        # it's easy to access any member of the parsed message
        #print secret.s.xs[0]  # first decoding basis coefficient of the secret

    A = []
    B = []

    instance = parse_instance(challenge, instance_path)
    if instance is None:
        print "Could not parse instance."
        sys.exit(-1)
    # print instance # print all contents
    # it is easy to access any member of the parsed message
    # each instance has many samples (print something from the first sample)
    # each sample is a pair (a,b=a*s+e) (print something in the ring element 'a')
    # each ring element is represented as a list of coefficients with resepct to
    #  the decoding basis. (print the first coefficient)
    # print instance.samples[0].a.xs[0]

    # print "See source code in ChallInstParser.py for examples of how to access specific message members."
    for i in range(challenge.dparams.numSamples):
        a = adjust_index(instance.samples[i].a.xs[0:])
        b = adjust_index(instance.samples[i].b.xs[0:])
        A.append(a)
        B.append(b)

    challID = challenge.challengeID
    n = challenge.dparams.m / 2
    q = challenge.dparams.q
    svar = challenge.dparams.svar
    numSample = challenge.dparams.numSamples
    return [challID,A, B, s, n, q, svar, numSample]
