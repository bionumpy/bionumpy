import sys
import re
import logging
logging.basicConfig(level=logging.INFO)

regex = r"cite\([a-z,0-9,A_Z, ]+\)"

with open(sys.argv[1]) as f:

    text = f.read()

    matches = re.finditer(regex, text , re.MULTILINE)

    for matchNum, match in enumerate(matches, start=1):

        logging.info("Match {matchNum} was found at {start}-{end}: {match}".format(matchNum=matchNum, start=match.start(),
                                                                            end=match.end(), match=match.group()))

        m = match.group()
        replacement = m.replace("cite(", ":cite:`").replace(")", "`")
        print("REplacement: " , replacement)
        text = text.replace(m, replacement)


print(text)