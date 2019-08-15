#!/usr/bin/python

# Parses RLWE secret files:
# takes the path to a .secret file and prints the contents.
# See the code for examples of how to access message contents.
#
# You will need the "protobuf" python package to run this code, which you can
# install with "easy_install protobuf"

# You can find more detailed information about parsing messages in python here:
# https://developers.google.com/protocol-buffers/docs/reference/python-generated
import sys
import Challenges_pb2

def parse_secret(secret_path):
  if secret_path is None:
    return
  secret = None
  with open(secret_path) as f:
    secret = Challenges_pb2.Secret()
    secret.ParseFromString(f.read())
  return secret

if __name__ == "__main__":

  if len(sys.argv) != 2:
    print "Usage:", sys.argv[0], "path/to/.secret"
    sys.exit(-1)

  secret = parse_secret(sys.argv[1])

  if secret is None:
    print "Could not parse secret."
    sys.exit(-1)
  else:
    print secret # print all contents
    # it's easy to access any member of the parsed message
    print secret.s.xs[0] # first decoding basis coefficient of the secret

  print "See source code in SecretParser.py for examples of how to access specific message members."

  sys.exit(0)