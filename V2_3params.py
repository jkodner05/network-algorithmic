# Extensions and parameter descriptions for V2 cf
# Niyogi & Berwick (1996) A langauge learning model for finite parameter spaces.
# Gibson & Wexler (1994) Triggers.

# 3 Parameter Model Gibson & Wexler 1994
paramdict = {}
#sfinal, cfinal, -v2
paramdict[1] = (1,1,0)
#sfinal, cfinal, +v2
paramdict[2] = (1,1,1)

#sfinal, cfirst, -v2
paramdict[3] = (1,0,0)
#sfinal, cfirst, +v2
paramdict[4] = (1,0,1)

#sfirst, cfinal, -v2
paramdict[5] = (0,1,0)
#sfirst, cfinal, +v2
paramdict[6] = (0,1,1)

#sfirst, cfirst, -v2
paramdict[7] = (0,0,0)
#sfirst, cfirst, +v2
paramdict[8] = (0,0,1)

sentencedict = {}
sentencedict[1] = set(["V S", "V O S", "V O1 O2 S", \
                   "Aux V S", "Aux V O S", "Aux V O1 O2 S", \
                   "Adv V S", "Adv V O S", "Adv V O1 O2 S", \
                   "Adv Aux V S", "Adv Aux V O S", "Adv Aux V O1 O2 S"])
sentencedict[2] = set(["S V", "S V O", "O V S", "S V O1 O2", "O1 V O2 S", "O2 V O1 S", \
                       "S Aux V", "S Aux V O", "O Aux V S", "S Aux V O1 O2", "O1 Aux V O2 S", "O2 Aux V O1 S", \
                       "Adv V S", "Adv V O S", "Adv V O1 O2 S", \
                       "Adv Aux V S", "Adv Aux V O S", "Adv Aux V O1 O2 S"])

sentencedict[3] = set(["V S", "O V S", "O2 O1 V S", \
                       "V Aux S", "O V Aux S", "O2 O1 V Aux S", \
                       "Adv V S", "Adv O V S", "Adv O2 O1 V S", \
                       "Adv V Aux S", "Adv O V Aux S", "Adv O2 O1 V Aux S"])
sentencedict[4] = set(["S V", "O V S", "S V O", "S V O2 O1", "O1 V O2 S", "O2 V O1 S", \
                       "S Aux V", "S Aux O V", "O Aux V S", "S Aux O2 O1 V", "O1 Aux O2 V S", "O2 Aux O1 V S", \
                       "Adv V S", "Adv V O S", "Adv V O2 O1 S", \
                       "Adv Aux V S", "Adv Aux O V S", "Adv Aux O2 O1 V S"])

sentencedict[5] = set(["S V", "S V O", "S V O1 O2", \
                       "S Aux V", "S Aux V O", "S Aux V O1 O2", \
                       "Adv S V", "Adv S V O", "Adv S V O1 O2", \
                       "Adv S Aux V", "Adv S Aux V O", "Adv S Aux V O1 O2"])
sentencedict[6] = set(["S V", "S V O", "O V S", "S V O1 O2", "O1 V S O2", "O2 V S O1", \
                   "S Aux V", "S Aux V O", "O Aux S V", "S Aux V O1 O2", "O1 Aux S V O2", "O2 Aux S V O1",\
                   "Adv V S", "Adv V S O", "Adv V S O1 O2", \
                   "Adv Aux S V", "Adv Aux S V O", "Adv Aux S V O1 O2"])

sentencedict[7] = set(["S V", "S O V", "S O2 O1 V", \
                       "S V Aux", "S O V Aux", "S O2 O1 V Aux", "Adv S V", \
                       "Adv S O V", "Adv S O2 O1 V", "Adv S V Aux", \
                       "Adv S O V Aux", "Adv S O2 O1 V Aux"])
sentencedict[8] = set(["S V", "S V O", "O S V", "S V O2 O1", "O1 V S O2", "O2 V S O1", \
                       "S Aux V", "S Aux O V", "O Aux S V", \
                       "S Aux O2 O1 V", "O1 Aux S O2 V", "O2 Aux S O1 V", \
                       "Adv V S", "Adv V S O", "Adv V S O2 O1", \
                       "Adv Aux S V", "Adv Aux S O V", "Adv Aux S O2 O1 V"])
