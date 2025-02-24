from cbor2 import load as load_data
from pprint import pprint
import sys

path = sys.argv[-1]
data = load_data(open(path, "rb"))
pprint(data)
