# CfD-QIS-IDL

Code to analyze and extract characteristic parameters (CG, RN, Offset, Quanta Exposure) of a Quanta Image Sensor (QIS). The method called Photon Counting Histogram (PCH) utilizes the photon-resolving capability of QIS. CG is calculated from the average distance between discrete peaks. RN, Offset, Quanta Exposure are extracted by curvefitting the QIS's histogram to a statistical model. QIS is currently being tested at room temperature and soon will be tested at cryogenic environment to better match the condition of space where its application is intended.

Learn more about QIS concept, development progress and future application: https://www.mdpi.com/1424-8220/16/8/1260
