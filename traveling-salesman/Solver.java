import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;

/**
 * The class <code>Solver</code> is an implementation of a greedy algorithm to solve the knapsack problem.
 *
*/

public class Solver {
  
    
  /**
  * The main class
  */
  public static void main(String[] args) {
    try {
      Solver solver = new Solver(args);
    } catch (IOException e) {
        e.printStackTrace();
    }
  }
    
  /**
  * Read the instance, solve it, and print the solution in the standard output
  */
  public Solver(String[] args) throws IOException {
    String fileName = null;
        
    // get the temp file name
    for(String arg : args){
      if(arg.startsWith("-file=")){
        fileName = arg.substring(6);
      } 
    }
    if(fileName == null)
      return;
        
    // read the lines out of the file
    List<String> lines = new ArrayList<String>();

    BufferedReader input =  new BufferedReader(new FileReader(fileName));
    try {
      String line = null;
      while (( line = input.readLine()) != null){
        lines.add(line);
      }
    }
    finally {
      input.close();
    }
              
    // parse the data in the file
    String[] firstLine = lines.get(0).split("\\s+");
    int N = Integer.parseInt(firstLine[0]);
    double x[] = new double[N];
    double y[] = new double[N];

    for (int i = 0; i < N; i++){
      String line = lines.get(i+1);
      String[] parts = line.split("\\s+");
      x[i] = Double.parseDouble(parts[0]);
      y[i] = Double.parseDouble(parts[1]);
    } 
   
    // execute the tsp on the input data
    TspSolver tsp = new TspSolver(x, y, false, false, true, true);
    //System.out.println(tsp.annealingParams());

    // output the solution
    if (!tsp.feasible()) System.out.println("infeasible");
    System.out.println(tsp);
/*   
    System.out.println("319857.4058964497 0");
    String s1 = "1576 1577 288 1363 1578 289 1579 1580 1581 1582 1583 1584 1064 1585 290 1722 696 713 599 1199 39 601 1868 521 937 331 930 1659 987 1853 931 702 739 459 623 391 552 514 134 693 932 506 131 588 324 353 347 332 980 975 1869 1198 1251 712 600 714 151 406 359 1682 356 1681 1364 355 1680 1670 1457 1679 1678 1677 775 575 495 463 682 936 718 917 266 543 425 941 570 555 603 889 640 874 656 501 884 618 878 625 891 518 875 549 86 85 1487 5 40 66 74 1851 46 84 49 1865 1888 1735 1836 1776 1844 681 486 1223 1222 1221 1220 1219 1218 1157 1156 1289 1428 1429 1430 1431 1432 1433 1434 1435 1436 1437 1438 1439 1440 1441 1442 1443 1444 1445 1446 1447 1448 1449 1450 1451 1452 1453 1454 1455 257 900 604 261 517 342 652 641 590 894 646 232 527 583 871 580 632 765 536 708 869 404 357 484 627 861 574 767 149 594 291 218 858 615 528 539 147 500 231 619 845 694 371 710 727 202 596 192 172 578 843 725 537 736 928 944 522 128 898 488 490 379 321 1845 981 1662 120 1660 988 1854 121 88 938 994 1863 703 740 460 1775 1817 1816 692 325 1674 982 329 532 497 345 926 351 132 326 589 507 933 690 550 135 457 621 515 737 700 939 392 993 990 989 986 983 1774 1772 1 991 1855 505 984 1852 1655 978 1588 1587 1197 1850 1040 1320 37 1059 25 1686 81 1060 38 1321 499 1586 1041 23 82 1687 107 29 1367 67 1872 31 918 1703 1799 1800 426 902 914 341 1669 1668 935 1798 424 423 264 1485 265 367 366 339 340 602 657 1473 1474 1475 1476 1477 1478 1479 1480 1481 1482 260 1483 890 676 1672 1871 30 915 903 233 28 1366 72 73 1685 52 24 1058 36 1319 43 1849 1362 1196 1292 95 977 976 1293 1195 1723 1361 1848 1045 529 1323 593 1057 56 55 1838 483 1684 1792 54 764 581 872 895 591 643 653 343 516 262 655 1864 642 1870 526 1365 1876 535 1837 1683 1791 668 667 534 1875 1671 1484 1759 1758 1488 1757 1756 1755 1754 1753 1752 1748 1472 1471 1470 1469 1468 1467 1466 1465 1464 1463 409 1793 1290 1159 1255 1019 1215 1462 1461 1523 1522 1521 307 1414 1413 1412 1411 1410 1409 1408 1407 1406 1405 1404 660 631 770 314 1598 1641 1640 1056 1639 1638 1637 963 964 965 162 966 419 386 1046 318 967 407 251 1259 1511 1644 1047 1645 1646 1600 316 1589 1647 1648 372 541 1599 547 393 731 1653 19 1617 1616 1416 436 1415 1628 1667 1458 375 1728 1700 1020 793 385 546 548 1794 1834 1783 1732 881 50 1733 1744 1743 485 1742 1741 1291 1885 1839 1835 1784 1734 1745 1867 1750 879 611 612 880 22 1882 27 893 892 885 886 883 882 876 877 21 554 4 553 1749 3 1866 0 1782 398 887 888 6 1881 26 1708 1709 1710 368 1711 1712 1713 1714 1715 1368 369 1716 1717 1718 1719 1720 1721 1360 1042 1322 1061 1688 70 69 647 896 1369 68 1873 32 34 1801 1841 1860 1704 1824 1705 1861 1842 613 1802 35 605 582 525 234 579 919 763 709 405 654 358 677 482 904 916 669 573 766 901 150 592 217 614 897 873 648 146 498 230 496 870 862 633 370 695 711 189 597 859 846 576 173 533 330 346 927 352 130 587 327 934 691 133 551 899 940 701 738 622 458 390 322 292 380 530 494 556 129 733 925 513 504 489 487 628 523 929 540 538 844 735 201 726 620 724 175 1194 1193 1192 1043 1191 1062 1190 1189 1188 1187 1186 170 1185 169 1184 1183 278 683 1239 1240 1241 1242 1243 1244 1245 1246 1063 1247 1044 1248 1249 191 1250 598 860 847 577 174 1843 979 985 1661 828 992 2 1773 48 827 190 671 672 670 443 1689 442 572 571 1884 1370 920 188 1238 1237 1236 1235 1234 1180 1847 1707 1883 1039 1038 1037 1036 144 1035 542 569 833 1880 832 835 119 834 1747 831 143 1034 1033 1032 1031 1030 1029 1028 1027 1026 1025 1024 145 1016 1015 389 1014 1731 1771 1460 1630 1770 1676 1419 1769 1665 1615 1420 1666 1502 1768 1652 1592 1603 1767 1766 1765 388 1051 1602 836 1591 373 319 323 794 568 472 271 333 304 336 337 686 474 299 334 250 328 338 303 335 707 595 544 394 734 272 477 456 249 354 723 445 662 309 503 449 438 259 924 732 606 374 401 427 365 524 607 752 481 256 182 216 193 142 263 415 1786 1785 1761 108 1524 1818 461 511 1746 624 921 51 922 626 1751 512 462 1819 1820 1821 1822 1823 1846 1706 1862 1179 1182 1181 168 171 200 1178 1177 1176 1175 1174 1173 1172 167 1171 199 198 166 1170 1169 197 165 1168 1167 1166 1165 1164 1163 1162 1161 1287 1286 1285 1023 1217 1284 1283 1282 1281 1280 1279 1278 1277 1276 1275 1274 1273 1272 1675 1418 722 1417 784 1629 1459 112 1730 1702 1216 1022 1739 1740 1160 1253 1252 1021 181 255 215 480 751 1701 1729 923 258 437 448 502 308 661 444 384 1605 706 300 1658 1501 1657 1656 90 20 1651 1726 1590 721 783 1601 705 1877 1777 1831 1048 968 1260 1815 1122 1536 782 704 1233 1830 719 720 252 254 1654 476 1261 196 969 1049 1262 1263 1264 1265 194 1266 1267 1268 1269 1270 1271";
    String  s2 = " 1317 1316 1315 1314 1313 1312 213 1311 1310 1309 1308 1050 970 214 471 344 148 387 408 479 520 420 315 317 127 212 195 273 531 475 320 253 447 161 473 279 226 610 717 220 187 470 637 634 363 685 617 312 311 305 645 639 154 680 452 277 378 140 431 675 455 689 564 301 559 562 772 269 468 210 248 811 754 825 294 584 276 821 774 268 422 759 867 943 237 863 853 383 280 757 817 848 838 780 428 699 786 842 852 746 829 857 778 756 204 814 636 805 792 810 789 178 185 567 466 493 1205 1206 1207 1208 1209 1210 1211 1212 1213 1214 180 560 13 1517 14 1595 1401 1307 16 11 1608 1003 1808 432 211 826 823 1494 12 179 1373 822 281 1542 1334 1541 1540 1539 913 1204 1333 1332 945 1538 1537 1696 1859 1001 912 1811 1790 1071 911 910 905 909 906 1070 1789 1810 907 1000 1858 1695 908 208 441 139 510 797 362 800 223 122 803 350 298 435 123 118 1694 1693 999 138 207 998 807 403 239 411 741 153 1692 1328 1329 1330 1331 949 1604 1825 1203 116 651 1202 788 184 177 465 492 809 566 297 125 434 802 799 349 222 509 361 796 206 440 137 1302 1788 240 1376 1069 1068 1301 997 1327 948 1384 1383 1382 1381 1380 1379 1378 1377 1067 241 1375 742 1066 63 1300 996 1857 1691 1326 98 99 1809 947 1226 1227 1201 1295 1228 1296 777 113 698 851 1388 1387 1386 1385 856 743 745 100 242 1389 1390 1391 1392 243 1393 247 1394 816 382 1371 1492 864 110 1515 275 1399 585 1398 246 1397 245 1396 1395 236 1807 1806 1805 1803 114 101 1297 64 1298 104 102 115 1804 111 105 1543 106 1372 716 91 1493 819 820 1516 33 1400 1306 15 10 1607 1002 44 1525 1072 1620 1635 45 61 1697 89 364 1634 1619 87 60 1699 92 1829 1232 609 225 160 1535 446 159 608 186 469 1698 684 1528 1527 1526 1007 1006 1005 1004 1609 1614 1305 1118 17 18 1495 1496 688 377 1499 1738 75 76 1613 53 1078 1625 451 1013 1621 954 1636 955 956 957 958 959 960 961 126 962 730 313 769 630 659 1510 1427 1258 1814 729 1121 1534 1337 1324 1231 1828 94 953 59 1531 1624 1633 1077 1012 83 1612 97 96 1737 1498 9 65 1011 1076 1632 1623 1530 58 952 93 1827 1230 219 42 1336 1533 1120 761 1813 1257 1426 1509 1124 1796 478 762 1763 1832 519 1125 1797 1833 1764 1673 47 1055 1650 1597 1879 1403 1887 1725 1594 1643 1500 1503 1618 1664 1421 1424 545 1760 1507 1627 1520 1519 1626 1506 80 1606 1423 1663 1778 1642 1593 1724 1886 1402 1878 1596 1649 1054 750 749 1762 666 665 1795 1123 1508 1425 1256 1812 748 664 1119 1532 1335 41 1229 1826 79 951 57 1529 1622 1631 1075 1010 679 1611 78 77 454 270 687 674 430 376 747 306 310 155 141 760 768 678 629 450 302 781 812 663 558 638 728 561 771 658 467 209 644 616 753 824 293 586 158 563 274 865 773 267 421 854 866 758 283 235 244 429 849 381 818 715 815 837 830 779 697 841 785 744 776 850 791 942 755 813 650 804 203 855 176 787 808 565 183 296 124 464 433 801 798 348 221 360 795 508 136 205 439 402 491 806 635 152 238 410 1374 1065 62 1787 1299 995 1856 1690 1325 117 295 946 649 1200 1294 71 790 103 839 840 282 1544 1545 1489 284 1546 1512 1115 1303 971 974 973 972 8 1304 1116 1513 7 557 1490 1491 868 109 1514 1117 453 673 1497 1736 1779 1780 1610 1009 1074 1008 1073 1079 1080 950 1081 1082 1083 224 227 416 412 395 156 1085 1084 285 399 1086 1087 1088 1089 1090 1091 400 286 1093 1092 157 396 413 417 228 1126 163 1127 164 1053 1128 1129 1130 1131 1132 1133 1134 1135 1136 1137 1138 1139 1140 1141 1142 1143 1144 1145 1146 1147 1148 1149 1150 1151 1152 1018 1153 1154 1155 1318 1358 1357 1356 1355 1354 1353 1352 1351 1350 1349 1348 1347 1346 1345 1344 1343 1342 1341 1340 1339 1338 229 418 414 397 1052 1094 1095 1096 1097 1098 1099 1100 1101 1102 1103 1104 1105 1547 1504 1548 1422 1549 1550 1505 1551 1518 1552 1553 1106 1107 1108 1109 1110 1111 1112 1113 1114 1017 1359 1554 1254 1158 1288 1727 1781 1874 1225 1840 1555 1556 1557 1224 1558 1559 1560 1561 1562 1563 1564 1565 1566 1567 1568 1569 1570 1486 287 1571 1572 1573 1574 1575 1456";
    System.out.println(s1 + s2); 
    */   
  }
  
}