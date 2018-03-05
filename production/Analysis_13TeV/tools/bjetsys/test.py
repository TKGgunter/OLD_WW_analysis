# coding: utf-8

import BTagCalibrationStandalone as BT
print "Intial test"
strings = [    "0, comb, central, 0, 0, 1, 0, 1, 0, 999, \"2*x\" \n",
    "0, comb, central, 0, 0, 1, 1, 2, 0, 999, \"2*x\" \n",
    "0, comb, central, 0, 1, 2, 0, 1, 0, 999, \"-2*x\" \n",
    "0, comb, central, 0, 1, 2, 1, 2, 0, 999, \"-2*x\" \n",
    "3, comb, central, 0, 0, 1, 0, 1, 2, 3, \"2*x\" \n",
    "3, comb, central, 0, -1, 0, 0, 1, 2, 3, \"-2*x\" \n",
    "1, test, central, 1, -2, 2, 50, 500, 0, 999, \"2*x\" \n",
    "1, test, up,      1, -2, 2, 50, 500, 0, 999, \"2.1*x\" \n",
    "1, test, down,    1, -2, 2, 50, 500, 0, 999, \"1.9*x\" \n"]

strings = u"".join(strings)

calibration = BT.BTagCalibration("csv")
calibration.readCSV(strings)

reader = BT.BTagCalibrationReader(BT.OperatingPoint.OP_MEDIUM, u"central", [u"up", u"down"])
reader.load(calibration, BT.JetFlavor.FLAV_C, u"test")

print "We expect 200, and we get ....", reader.eval_auto_bounds(u"central", BT.JetFlavor.FLAV_C, 0.5, 100., 0)
print BT.bjetScale(u"central",  [.23], [4], [0.5], [20], reader)





print "loading data test"
f = open("data/CSVv2_Moriond17_B_H.csv")
strings = f.read()
f.close()


print "readCSV"
calibration = BT.BTagCalibration("csv")
calibration.readCSV(strings)


print "loading calibration"
reader = BT.BTagCalibrationReader(BT.OperatingPoint.OP_MEDIUM, u"central", [u"up", u"down"])
reader.load(calibration, BT.JetFlavor.FLAV_B, u"comb")


print reader.eval_auto_bounds(u"central", BT.JetFlavor.FLAV_B, 0.5, 100., 0)
print BT.bjetScale(u"central",  [0.8485]*9, [5]*9, [0.5]*9, [21 + i for i in range(9)], reader)
print BT.bjetScale(u"up",       [0.8485]*9, [5]*9, [0.5]*9, [21 + i for i in range(9)], reader)
print BT.bjetScale(u"down",     [0.8485]*9, [5]*9, [0.5]*9, [21 + i for i in range(9)], reader)

