package tsa;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.app.util.TreeFile;
import beast.app.util.XMLFile;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Logger;
import beast.core.Operator;
import beast.core.Runnable;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.datatype.Binary;
import beast.evolution.likelihood.AncestralStateTreeLikelihood;
import beast.evolution.likelihood.ThreadedTreeSetLikelihood;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.operators.DeltaExchangeOperator;
import beast.evolution.operators.IntRandomWalkOperator;
import beast.evolution.operators.ScaleOperator;
import beast.evolution.operators.SwapOperator;
import beast.evolution.operators.UniformOperator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.gss.NS;
import beast.math.distributions.Gamma;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.Prior;
import beast.util.TreeParser;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import beast.util.XMLProducer;
import tsa.correlatedcharacters.polycharacter.CompoundAlignment;
import tsa.correlatedcharacters.polycharacter.CorrelatedSubstitutionModel;
import tsa.parameterclone.helpers.RescaledDirichlet;
import tsa.parameterclone.selector.Selector;

@Description("Perform correlated character analysis on D-PLACE data")
public class TSAModelSelector extends Runnable {
	public Input<Integer> particleCountInput = new Input<>("particleCount", "number of particles (default 1)", 1);
	public Input<Integer> subChainLengthInput = new Input<>("subChainLength",
			"number of MCMC samples for each epoch (default 100)", 100);
	public Input<File> dataFileInput = new Input<>("file", "CSV file exported from DPLACE data",
			new File("dplace/dplace-societies-2017-02-21.csv"));
	public Input<String> thresholdV1Input = new Input<>("v1", "Comma separated list of values for first DPLACE feature "
			+ "that are interpreted as 1. "
			+ "For example \"1,2,4\" means when the feature's value is 1 or 2 or 4 it becomes a 1, and otherwise a 0 "
			+ "(but \"NA\" is interpreted as missing data)"
			, Validate.REQUIRED);
	public Input<String> thresholdV2Input = new Input<>("v2", "As thresholdV1 but for second DPLACE feature",
			Validate.REQUIRED);
	public Input<Double> epsilonInput = new Input<>("epsilon",
			"stopping criterion: smallest change in ML estimate to accept", 1e-8);
	public Input<OutFile> outputInput = new Input<>("output", "where to save the dot file with results",
			new OutFile("/tmp/x.dot"));
	public Input<XMLFile> xmlFileInput = new Input<>("xml", "XML file containing NS analysis, so you can set your own priors. "
			+ "The IntegerParameter with id=\"indices\" will be replaced by one of the 16 possible combinations.");

	public Input<TreeFile> treeSetFileInput = new Input<>("treeFile", "file containing trees used for inference");
	
	public Input<Integer> threadCountInput = new Input<>("threads", "Number of threads to use (default 1)", 1);
	public Input<Boolean> useRelaxedClockInput = new Input<>("relaxed", "if true, a relaxed clock is used, otherwise a strict clock is used",
			false);

	public Input<String> codeInput = new Input<>("code", "column name containing the code for taxa in the tree, e.g  or 'Society id' or 'ISO code'", "Society id");

	public Input<XMLFile> exportXmlFileInput = new Input<>("export", "If specified, XML file is exported containing NS analysis that is run.");

	
	@Override
	public void initAndValidate() {
	}

	String[] indices;
	double[] ml;
	double[] df;
    public static ExecutorService exec;
	CountDownLatch countDown;
	DecimalFormat format = new DecimalFormat("##0.0");
	int [] siteCount;

	final static String newick = "((((((((((xd428,xd425),((xd432,xd433),xd430)),xd431),xd378),xd441),xd440),((xd492,xd498),(xd467,xd468))),(xd374,(((xd382,xd507),xd460),xd466))),(((((((xd442,xd453),(xd443,xd452)),(((xd454,((((xd424,xd427),(xd416,xd426)),((xd417,xd419),xd415)),xd423)),(xd420,(xd421,xd422))),((xd445,xd439),(xd457,xd448)))),((xd590,(((((((((xd515,xd512),(xd494,xd510)),(xd513,xd514)),((xd517,xd520),xd511)),(xd487,((xd495,(xd504,xd497)),(xd527,xd496)))),xd519),((((((xd463,xd465),((xd481,((xd503,xd461),xd462)),xd464)),xd516),xd518),xd588),((xd587,((xd586,(xd593,xd594)),xd592)),(xd595,xd591)))),xd589),((xd455,xd456),(xd437,(xd446,xd444))))),(((((((xd476,(xd486,xd458)),xd475),((xd478,(xd459,xd484)),xd305)),(xd307,xd479)),xd383),(((((xd477,xd325),xd327),(xd328,xd483)),xd334),((((xd482,xd347),xd346),(((xd326,xd332),xd330),xd345)),xd329))),(((((xd525,(xd523,xd493)),xd522),((xd526,(xd508,xd501)),xd521)),(xd509,xd499)),(((xd489,xd506),((xd505,xd488),(xd490,(xd502,xd491)))),xd500))))),xd524),(xd531,xd535)),(((((((((((((((((xd1390,xd1395),xd1391),xd1393),(xd1374,xd1396)),((xd1371,xd1382),((xd1403,(xd1392,xd1388)),xd1401))),xd919),(xd1375,xd1399)),((xd922,xd1384),(((xd1385,xd1386),xd1383),xd1387))),((((((((((((((xd920,xd921),xd1408),xd1394),(xd1362,xd918)),((xd1363,xd1416),(xd1402,(xd1354,xd1414)))),xd1406),xd1357),xd1353),(xd1397,xd1405)),(((((xd1417,xd1410),xd1409),(xd1415,xd1412)),xd1411),(xd1407,xd1413))),(((((((xd916,xd1342),(xd1329,(xd1346,xd1326))),xd1340),(xd1400,xd1398)),(((((xd1337,xd1341),xd1338),xd1339),xd1350),xd1359)),xd1332),((xd1358,(xd1361,xd1356)),xd1360))),(((((((xd1328,xd1377),xd1379),xd1327),((xd1309,xd1314),xd1321)),xd1347),(xd1370,xd917)),((xd1369,xd1365),(xd1351,xd1335)))),((((xd1372,xd1364),xd1373),xd1368),xd1345)),(((((xd1336,((xd1323,(xd1330,xd1331)),xd1333)),((xd1352,xd1343),xd1348)),((((xd1380,xd1381),xd1389),xd1404),xd1355)),xd1320),((xd1376,xd1367),xd1366)))),xd1378),(((((((((xd1291,xd1287),xd1294),(((xd1286,xd1293),((xd1308,xd1301),xd1288)),xd1292)),((((((xd1200,(xd1158,(((xd1157,xd1154),xd1155),(xd1156,xd994)))),(((xd989,xd1161),xd992),xd1160)),(xd1122,((xd993,xd1162),xd1164))),((xd1190,((xd1198,(xd1184,((xd1199,xd1193),((xd1194,(xd1152,(xd1197,xd1186))),(xd1195,xd991))))),(xd1188,((xd1189,xd1192),xd1187)))),xd1126)),(xd1222,(((xd1183,(xd1182,((xd1181,xd1185),(xd1180,xd1201)))),(xd1168,((xd1169,(xd1178,xd1171)),((xd995,((xd1175,xd1179),((xd988,xd1176),xd1177))),(xd1172,xd1170))))),(xd1166,(xd1165,xd1167))))),(((((xd1132,((xd1439,xd1440),xd1131)),xd1133),xd1118),(xd1130,xd1128)),xd1268))),((xd1107,xd1112),xd1121)),((((xd1094,((((xd1269,(xd1437,(xd1137,xd1271))),((xd1281,xd1273),xd1272)),(((xd1163,xd1202),(xd1203,(xd1204,xd1274))),xd1426)),xd1134)),((xd1111,xd1110),(xd1433,xd1432))),xd1109),(xd1127,xd1289))),xd1316),((((xd1270,(((xd1430,xd1277),xd1429),xd1261)),xd1284),(xd1216,((xd1285,xd1282),xd1260))),xd1280)),((((xd1254,(xd1245,xd1244)),xd1215),xd1095),xd1093))),(((((((((xd1108,xd1135),((((xd1283,xd1262),xd1278),xd1276),(xd1263,xd1264))),((((xd1209,xd1240),xd1239),(xd1242,xd1234)),((xd1078,xd1150),(((xd1435,xd1149),(xd1141,xd1151)),xd1153)))),(xd1081,((((xd1129,xd1436),xd1123),xd1124),(xd1140,xd1139)))),((((((((((((xd1049,(xd1050,xd1053)),xd1213),((xd1051,xd1055),(xd1052,xd1054))),(xd1004,(xd1003,xd1006))),xd1057),xd1059),(((xd1235,xd1241),xd1233),xd1238)),((xd1058,(xd1253,xd1248)),((xd1048,xd1063),((xd1218,xd1065),(((xd1062,xd1047),xd1005),xd1011))))),(((xd1212,(xd1211,xd1217)),xd1224),xd1243)),(xd1228,xd1208)),(xd1097,xd1101)),(((xd1088,((((xd1106,(xd997,xd1085)),xd1080),(((((((xd1087,xd1086),xd1441),(xd1427,xd1079)),(xd1073,xd1071)),xd1072),((xd1075,(xd1074,xd1428)),xd1076)),((((xd1146,(xd1419,xd1148)),(xd1145,xd1144)),(xd1174,(xd1001,xd1147))),(xd1207,(xd1143,xd1142))))),(((((xd1083,xd1424),xd1090),xd1105),(xd1084,xd1070)),((xd998,xd1104),xd1077)))),(xd1205,xd1125)),((((xd1136,xd1421),xd1420),(xd983,xd1120)),xd1089)))),(((((xd1028,xd622),(xd1027,(xd1014,(xd1064,xd1000)))),(((xd1029,(xd1036,(xd1020,xd1021))),((xd1015,xd1016),xd1018)),((((xd1030,(xd1039,(xd1061,xd1017))),((xd1040,xd1041),xd1031)),(((xd1023,xd1060),(xd1038,xd1024)),(xd1022,xd1019))),xd1025))),(xd999,xd1067)),(((((((xd1045,xd1044),(xd1046,(xd1032,(xd1034,((xd1033,xd1009),xd1007))))),((((xd1035,xd1069),(xd1056,xd1043)),(xd1206,xd1226)),(((((xd1066,(xd1013,xd1010)),xd1042),xd1037),xd1026),(xd1008,xd1012)))),((((xd1423,xd1092),(xd1103,xd1422)),(((xd1096,xd1425),xd1099),(xd1100,(xd1434,xd1098)))),(((xd1266,xd1275),((xd1267,(xd1259,xd1265)),xd1279)),xd1219))),xd1102),(xd1431,xd1082)),(xd1068,xd1438)))),xd1091),((((((((xd1210,(xd1220,xd1227)),(xd1237,xd1232)),(xd1236,(xd1230,xd1231))),((xd1223,xd1214),xd1225)),xd1258),xd1246),xd1196),xd1256)),(((((xd1290,xd1299),xd1295),((xd1304,xd1296),xd1298)),(((((xd1319,xd1318),xd1315),(xd1310,xd1322)),(((xd1311,xd1313),xd1317),xd1324)),xd1297)),(xd1300,((xd1303,(xd1302,xd1306)),(xd1307,xd1305)))))),(((((xd1249,xd1247),(xd1252,xd1251)),xd1221),((xd1250,xd985),(xd1119,((xd1113,xd1114),xd1138)))),((xd1116,((xd1418,xd1117),xd1115)),((xd1312,xd1325),xd1257)))),((xd631,xd633),xd621)),((((((((xd577,xd580),xd654),((xd584,xd578),xd585)),((((((((xd537,xd541),(xd1349,(xd536,xd538))),xd542),((((xd547,xd563),((((xd539,xd529),xd533),((xd1334,xd540),xd545)),xd532)),xd530),xd543)),(((((xd570,xd548),(xd571,xd557)),xd567),(xd569,xd546)),(xd559,xd573))),(xd528,xd534)),(((((xd607,((xd605,xd608),xd602)),((xd576,xd597),xd582)),(xd598,xd601)),(((((((xd664,xd697),xd683),(xd663,xd662)),(xd692,xd690)),(((xd667,xd671),xd685),(((xd603,xd657),(xd670,(xd600,xd596))),xd669))),xd665),xd604)),xd575)),xd583)),((((((((xd549,(xd558,xd561)),xd609),(xd560,xd579)),(xd616,(xd599,xd613))),(xd554,xd581)),xd628),xd551),((((xd568,(xd552,((((xd553,xd562),(xd555,xd566)),xd556),xd565))),xd572),(((xd618,xd623),(xd550,xd564)),xd544)),((xd632,xd620),xd619)))),(((xd612,xd610),((((xd611,xd615),xd614),xd574),xd606)),xd636)),((((xd624,(xd625,xd915)),((xd637,(xd627,xd626)),(xd629,xd630))),xd646),xd634)),((xd686,xd901),(((((((xd672,xd677),xd673),((xd679,xd674),xd676)),xd675),((((((xd912,xd910),xd682),(xd905,xd678)),xd681),xd680),(xd913,xd914))),xd668),(xd903,xd902))))),((xd635,xd617),(((xd650,(xd642,xd640)),xd648),xd638))),((((xd764,(((((((((((((xd803,xd805),(xd812,xd810)),(xd806,xd818)),((((xd848,xd851),xd857),xd846),(xd852,(xd799,xd786)))),(((xd847,xd859),(xd856,(xd860,(xd855,xd844)))),xd849)),(((((xd871,(xd865,xd863)),(((xd864,(xd866,xd869)),xd868),xd862)),xd873),(((xd850,xd858),xd853),xd875)),(((((xd829,xd838),(((((xd830,(xd833,xd825)),(xd831,xd835)),xd840),xd836),(xd827,(xd832,xd834)))),xd826),xd828),((((xd867,xd874),xd870),(((((xd887,(xd877,xd881)),(xd882,xd879)),((xd885,(xd880,xd876)),xd883)),((((((xd895,((xd890,xd893),(xd897,xd894))),(xd891,xd889)),(xd896,xd892)),xd898),xd884),xd886)),(xd888,xd878))),xd872)))),(((xd842,xd861),xd845),xd837)),xd820),((((xd771,xd763),xd769),(xd770,(xd766,xd765))),xd768)),((((((((xd907,xd906),xd908),(((xd748,xd742),xd743),xd739)),(xd904,xd737)),(((xd741,xd735),((xd900,xd747),xd745)),(xd744,xd738))),xd740),(xd839,xd824)),((((xd761,(xd750,((xd757,xd730),(xd755,xd723)))),(xd718,xd717)),(xd752,xd762)),((((xd753,xd756),xd728),((((xd754,xd1442),xd758),(xd759,xd767)),xd751)),((((((xd687,xd688),xd689),xd695),xd694),(xd928,xd693)),xd736))))),xd746),(xd749,(xd732,xd733))),xd734)),xd854),(((((((((xd698,(xd710,xd706)),(xd708,(xd712,xd699))),(((xd700,(xd702,xd703)),(xd701,xd707)),xd704)),((xd653,(xd658,xd656)),(xd660,xd655))),((xd696,xd659),xd711)),((xd729,xd645),(xd709,xd715))),(xd651,(((xd639,xd643),xd644),xd649))),xd713),(((((xd911,(xd666,xd684)),xd661),((xd724,((xd725,xd721),(xd722,xd727))),(xd909,(xd720,xd726)))),(xd691,(xd714,((xd716,xd719),xd705)))),((xd652,xd731),(xd647,xd641))))),(((((((((((((((((xd964,xd772),xd968),(xd974,xd972)),((xd777,xd940),(xd779,xd963))),((((xd956,xd773),xd953),xd971),((xd958,xd959),((xd965,(xd967,xd969)),xd960)))),(((((xd947,xd781),xd945),xd775),(xd946,xd955)),xd937)),(xd780,xd976)),((((xd966,xd961),(((xd977,xd978),xd975),(xd970,xd784))),((xd973,xd979),xd962)),(xd954,xd957))),(xd952,xd951)),xd948),(((xd776,xd944),xd939),xd782)),xd950),(((xd941,xd774),(xd811,xd814)),(xd794,xd822))),((xd938,xd789),xd778)),(((((xd804,(xd807,xd790)),(xd792,xd843)),(((xd802,xd841),xd817),xd943)),((((xd785,xd816),(xd809,xd793)),(((xd821,(xd787,xd796)),(xd808,(xd815,xd800))),xd823)),((xd801,xd797),xd760))),xd980)),((((xd791,xd813),(xd949,xd942)),xd798),((xd788,xd819),xd795))),xd783))))),(((((((((((((((xd217,(xd189,(xd208,xd222))),(xd218,xd243)),xd220),((xd224,xd225),xd221)),xd216),((xd215,xd274),xd186)),(xd214,xd223)),((xd226,xd241),(xd228,(xd229,xd227)))),xd340),((((((((xd201,xd187),(((xd202,xd204),xd203),(xd190,(xd206,xd242)))),(xd314,((xd313,xd309),xd315))),((((xd196,xd198),xd219),xd207),(((((xd211,xd240),xd212),xd213),xd200),xd209))),((((xd299,((xd333,xd311),xd336)),(xd300,xd301)),(xd303,xd302)),xd304)),(((xd210,xd306),((xd205,xd312),xd308)),(xd310,xd298))),(((((((((xd127,(xd154,(xd150,(xd5,xd152)))),xd153),xd151),((((((xd169,xd149),xd162),xd175),((xd168,xd130),(xd138,xd166))),(((xd164,xd161),xd163),(xd167,xd40))),(xd165,xd173))),(((((((xd105,(((xd92,((xd122,xd96),xd77)),xd121),(xd95,xd81))),((xd76,xd97),((xd99,(xd100,xd101)),xd102))),((xd98,xd106),((xd88,xd90),xd117))),((((((xd78,xd125),xd91),(xd65,xd75)),xd85),(((((xd111,xd126),xd110),(((((xd80,xd115),(xd87,xd118)),xd112),(xd113,xd119)),(xd124,xd116))),(xd79,xd107)),((((xd132,xd128),xd184),(xd133,xd131)),xd114))),((xd103,(xd109,xd108)),xd104))),(((((((((xd24,(xd11,xd12)),xd13),xd21),xd32),xd27),(xd25,xd17)),(((((xd15,xd16),(xd26,xd31)),xd14),xd69),xd29)),((((xd59,xd71),(xd60,xd55)),xd61),xd94)),(((xd57,xd925),xd129),xd135))),(((((((((xd34,xd64),xd67),xd63),(xd37,xd62)),(xd18,xd19)),((((xd86,(xd70,xd82)),xd89),xd84),(xd83,xd123))),(xd68,xd73)),(xd93,xd120)),(xd33,xd56))),((((xd923,xd924),((xd146,xd145),xd156)),((((xd58,xd137),xd182),xd134),xd136)),((((((xd159,xd158),xd148),(xd157,xd147)),(((((xd66,xd48),xd49),xd140),(xd141,xd144)),((xd139,(xd143,xd160)),xd142))),(((xd47,xd51),(xd46,xd50)),xd43)),((((((xd36,xd35),(xd52,xd53)),xd72),((xd10,(xd20,(xd22,xd30))),xd28)),(((xd39,xd54),xd44),xd23)),(((xd38,xd74),xd42),(xd41,xd45))))))),((xd317,xd320),xd195)),((((xd185,((xd172,xd170),xd179)),(xd171,xd183)),xd181),((xd176,xd174),xd178))),(xd194,xd193)),((xd177,xd180),(((xd331,xd335),xd316),xd319)))),(((((((xd926,xd352),((xd356,xd359),xd355)),(xd354,xd337)),((xd341,xd342),(xd343,xd349))),(((xd344,xd323),xd318),xd322)),(xd353,xd381)),((xd272,(((xd275,(((((xd282,xd284),xd279),(((xd273,xd288),xd276),((xd247,xd277),(xd255,xd285)))),(((xd283,xd289),xd245),xd286)),xd290)),(xd321,xd324)),(xd271,((xd270,xd246),(((xd287,xd291),(xd278,xd281)),(xd248,xd280)))))),(xd269,xd268))))),(((((xd471,((xd472,(xd485,xd474)),(xd294,xd473))),xd258),((xd251,xd469),((xd257,xd296),xd252))),((xd250,(xd254,xd256)),xd253)),(((xd249,xd238),(xd197,xd191)),xd236))),xd199),((((((((((xd392,xd411),(xd385,xd398)),xd409),((xd375,xd408),xd380)),xd376),(xd406,xd386)),((xd404,((((xd399,xd402),xd2),((xd401,xd412),xd414)),xd403)),((((((xd397,xd407),xd410),xd400),xd384),((xd395,xd405),xd387)),((xd389,xd391),(xd390,(xd413,xd388)))))),((xd394,(xd393,xd396)),xd370)),((((((xd350,xd351),xd365),xd366),(((xd927,xd360),xd379),(((xd361,xd362),xd363),xd364))),(xd358,(xd155,(xd339,xd357)))),(xd377,xd451))),((((((xd434,xd438),((xd435,(xd450,xd429)),(xd449,xd418))),(xd447,xd436)),(((xd367,xd338),(xd371,xd372)),xd348)),xd369),(xd368,xd373)))),(((((((((xd244,xd263),(xd297,xd259)),xd260),xd239),(xd262,xd293)),(((xd235,xd230),xd188),xd192)),(((xd264,xd295),xd261),xd267)),((xd292,xd265),((((xd237,xd234),xd231),xd233),xd232))),(xd266,(xd480,xd470)))),((((((xd4,xd3),xd930),(((xd932,xd7),xd929),xd933)),(xd934,((xd935,xd8),xd936))),(xd1,xd931)),(xd9,xd6))));";

	
	@Override
	public void run() throws Exception {
		long start = System.currentTimeMillis();
        Log.warning("Do not use BEAGLE with CorrelatedSubstitutionModel: setting java.only=true");
        System.setProperty("java.only", "true");

        indices = new String[16];
		indices[0] = "0 1 2 1 3 0 3 2";

		indices[1] = "0 1 2 1 3 0 4 2";
		indices[2] = "0 1 2 1 3 0 3 4";
		indices[4] = "0 1 2 4 3 0 3 2";
		indices[8] = "0 1 2 1 3 4 3 2";

		indices[3] = "0 1 2 1 3 0 4 5";
		indices[5] = "0 1 2 4 3 0 5 2";
		indices[6] = "0 1 2 3 4 0 4 5";
		indices[9] = "0 1 2 1 3 4 5 2";
		indices[10] = "0 1 2 1 3 4 3 5";
		indices[12] = "0 1 2 3 4 5 4 2";

		indices[7] = "0 1 2 3 4 0 5 6";
		indices[11] = "0 1 2 1 3 4 5 6";
		indices[13] = "0 1 2 3 4 5 6 2";
		indices[14] = "0 1 2 3 4 5 4 6";

		indices[15] = "0 1 2 3 4 5 6 7";

		ml = new double[16];
		df = new double[16];

		int threadCount = threadCountInput.get();
		if (threadCount == 1) {
			for (int i = 0; i < 16; i++) {
				Log.warning("Processing model " + i + ": " + indices[i]);
				NS NS = buildModel(indices[i]);
				NS.run();
				ml[i] = NS.getMarginalLikelihood();
				df[i] = NS.getStandardDeviation();
			}
		} else {
		    exec = Executors.newFixedThreadPool(threadCount);
			countDown = new CountDownLatch(16);
			// kick off the threads
			for (int i = 0; i < 16; i++) {
				CoreRunnable coreRunnable = new CoreRunnable(i);
				exec.execute(coreRunnable);
			}
			countDown.await();
		}

//		ml = new double[] { -711.8 , -702.5 , -709.6 , -703.3 , -699.6 , -696.2 , -705.0 , -709.8 , -704.9 , -720.0 , -702.9 , -706.5 , -707.1 , -702.1 , -698.7 , -711.5};
		reportStats();		
		
		PrintStream out = System.out;
		if (outputInput.get() != null) {
			out = new PrintStream(outputInput.get());
			Log.info("\nOutput written to " + outputInput.get().toPath() +"\n");
		}
		String dot = toDot(ml, df);
		out.print(dot);

		long end = System.currentTimeMillis();
		Log.warning("Total time spent: " + (end - start) / 1000 + " seconds");
		
		Log.warning("All done!");

	}

	private void reportStats() {
		Log.info("");
		Log.info("model  : indices        : marginal likelihood (standard deviation)");
		for (int i = 0; i < 16; i++) {
			String bin = Integer.toBinaryString(i);
			while (bin.length() < 4) {
				bin = "0" + bin;
			}
			Log.info((i < 10 ? " ":"") + i + " " + bin + ": " + indices[i] + " : " + format.format(ml[i]) + " (" + format.format(df[i]) + ")");
		}		
		Log.info("");
		Log.info("rate    : P(equal) P(different)");
		
		odds(1, "?0 -> ?1");
		odds(2, "0? -> 1?");
		odds(3, "?1 -> ?0");
		odds(4, "1? -> 0?");		
		
		Log.info("");
		Log.info("value: count");
		Log.info("00   : " + siteCount[0]);
		Log.info("01   : " + siteCount[1]);
		Log.info("10   : " + siteCount[2]);
		Log.info("11   : " + siteCount[3]);
	}

	private void odds(int rate, String label) {
		double max = ml[0];
		for (int j = 0; j < 16; j++) {
			max = Math.max(ml[j], max);
		}
		int [] isOn = new int[16];
		switch (rate) {
		case 1: isOn = new int[]{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};break;
		case 2: isOn = new int[]{0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1};break;
		case 3: isOn = new int[]{0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1};break;
		case 4: isOn = new int[]{0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};break;
		}
		double [] p = new double[2];
		for (int j = 0; j < 16; j++) {
			p[isOn[j]] += Math.exp(ml[j] - max);
		}
		double sum = (p[0] + p[1]) / 100.0;
		p[0] /= sum;
		p[1] /= sum;
		Log.info(label + ": " + 
				"      ".substring(format.format(p[0]).length()) + format.format(p[0]) + "%" + 
				"          ".substring(format.format(p[1]).length()) + format.format(p[1]) + "%");
	}




	class CoreRunnable implements java.lang.Runnable {
		int i;

		CoreRunnable(int i) {
			this.i = i;
		}

		@Override
		public void run() {
			try {
				Log.warning("Processing model " + i + ": " + indices[i]);
				NS NS = buildModel(indices[i]);
				NS.run();
				ml[i] = NS.getMarginalLikelihood();
				df[i] = NS.getStandardDeviation();
			} catch (Exception e) {
				e.printStackTrace();
			}
			countDown.countDown();
		}

	} // CoreRunnable

	private String toDot(double[] ml, double[] df) {
		String dot = "digraph {\n" + " graph [mindist=0.0, nodesep=0.25, ranksep=0.4]\n"
				+ "; node [fontsize=\"12\", style=\"solid\", color=\"#0000FF60\", width=0.95, height=0.75, fixedsize=\"true\",color=\"#FFFFFF\"];\n"
				+ " \n" + "\n" + "0 -> 1;\n" + "0 -> 2;\n" + "0 -> 4;\n" + "0 -> 8;\n" + "\n" + "1 -> 3;\n"
				+ "1 -> 5;\n" + "1 -> 9;\n" + "2 -> 3;\n" + "2 -> 6;\n" + "2 -> 10;\n" + "4 -> 5;\n" + "4 -> 6;\n"
				+ "4 -> 12;\n" + "8 -> 9;\n" + "8 -> 10;\n" + "8 -> 12;\n" + "\n" + "3 -> 7;\n" + "3 -> 11;\n"
				+ "5 -> 7;\n" + "5 -> 13;\n" + "6 -> 7;\n" + "6 -> 14;\n" + "9 -> 11;\n" + "9 -> 13;\n" + "10 -> 11;\n"
				+ "10 -> 14;\n" + "12 -> 13;\n" + "12 -> 14;\n" + "\n" + "7 -> 15;\n" + "11 -> 15;\n" + "13 -> 15;\n"
				+ "14 -> 15;\n" + "\n";
		for (int i = 0; i < 16; i++) {
			String bin = Integer.toBinaryString(i);
			while (bin.length() < 4) {
				bin = "0" + bin;
			}
			dot += i + " [label=\"" + bin + "\\n" + format.format(ml[i]) + " (" + format.format(df[i]) + ")\"]\n";
		}
		dot += "\n" + "\n" + "}\n";
		return dot;
	}
	
	private void toXML(NS o) throws IOException {
		if (exportXmlFileInput.get() != null) {
			Logger logger = new Logger();
			logger.setID("screenlog");
			logger.everyInput.setValue(1000, logger);
			logger.loggersInput.setValue(o.posteriorInput.get(), logger);
			o.loggersInput.setValue(logger, o);

			XMLProducer p = new XMLProducer();
			String xml = p.toXML(o);
	        FileWriter outfile = new FileWriter(exportXmlFileInput.get());
	        outfile.write(xml);
	        outfile.close();
	        
	        o.loggersInput.get().clear();
		}		
	}

	private NS buildModel(String indexValues) throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		if (xmlFileInput.get() != null && !xmlFileInput.get().equals("[[none]]")) {
			NS o = buildModelFromFile(indexValues);
			toXML(o);
			return o;
		}
		if (useRelaxedClockInput.get()) {
			NS o = buildModelRelaxed(indexValues);
			toXML(o);
			return o;
		} else {
			NS o = buildModelStrict(indexValues);
			toXML(o);
			return o;
		}
	}
	
	private NS buildModelFromFile(String indexValues) throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		XMLParser parser = new XMLParser();
		Object o = parser.parseFile(xmlFileInput.get());
		if (!(o instanceof NS)) {
			throw new IllegalArgumentException("Expected the run element to be of type NS");
		}
		NS ns = (NS) o;
		IntegerParameter indices = getIndicesObject(ns);
		if (indices == null) {
			throw new IllegalArgumentException("Expected element with id=\"indices\" in XML, but could not find any");
		}
		indices.valuesInput.setValue(indexValues, indices);
		indices.initAndValidate();
		return ns;
	}

	private IntegerParameter getIndicesObject(BEASTInterface o) {
		for (BEASTInterface bi : o.listActiveBEASTObjects()) {
			if (bi.getID().equals("indices")) {
				if (bi instanceof IntegerParameter) {
					return (IntegerParameter) bi;
				} else {
					throw new IllegalArgumentException("Element with id=\"indices\" is not an IntegerParameter");
				}
			} else {
				IntegerParameter result = getIndicesObject(bi);
				if (result != null) {
					return result;
				}
			}
		}
		return null;
	}

	private NS buildModelStrict(String indexValues) throws IOException {
		TreeParser tree = new TreeParser(null, newick, 0, false);

		Alignment data = processDataFile(tree.getTaxaNames());
		CompoundAlignment characters = new CompoundAlignment(data);

		// <run id="mcmc" spec="beast.gss.NS" chainLength="50000" preBurnin="5"
		// particleCount="1" subChainLength="100">

		// <state id="state" storeEvery="100">
		//
		// <stateNode spec="beast.core.parameter.RealParameter" id="parameters"
		// lower="0.05">
		// 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
		// </stateNode>
		RealParameter parameters = new RealParameter();
		parameters.initByName("value", "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0", "lower", 0.05, "upper", 8.0);
		parameters.setID("parameters");
		// <stateNode spec="beast.core.parameter.IntegerParameter" id="indices"
		// lower="0" upper="2">
		// 0 1 2 3 4 5 6 7
		// </stateNode>
		IntegerParameter indices = new IntegerParameter(indexValues);
		indices.setID("indices");
		// <stateNode spec="beast.core.parameter.IntegerParameter" id="sizes">
		// 1 1 1 1 1 1 1 1
		// </stateNode>
		IntegerParameter sizes = new IntegerParameter("1 1 1 1 1 1 1 1");
		sizes.setID("sizes");
		// <parameter id="frequencies" name="stateNode" lower="0" upper="1"
		// value="0.25 0.25 0.25 0.25"/>
		RealParameter frequencies = new RealParameter("0.25 0.25 0.25 0.25");
		frequencies.setID("frequencies");
		// <stateNode idref="clockRate.c:dplace"/>
		RealParameter clockRate = new RealParameter("0.5");
		clockRate.initByName("upper", 1.0, "value", "0.5");
		clockRate.setID("clockRate");
		// </state>

		State state = new State();
		state.initByName("stateNode", parameters,
				// "stateNode", indices,
				"stateNode", sizes,
				// "stateNode",frequencies,
				"stateNode", clockRate);

		// <operator id="DeltaExchanger" spec="DeltaExchangeOperator"
		// weight="1">
		// <parameter name="weightvector" idref="sizes" />
		// <parameter name="parameter" idref="parameters" />
		// </operator>

		DeltaExchangeOperator deltaExchanger = new DeltaExchangeOperator();
		deltaExchanger.initByName("weightvector", sizes, "parameter", parameters, "autoOptimize", false, "weight", 1.0);
		//
		// <operator id="clockRateScaler.t:dna" spec="ScaleOperator"
		// scaleFactor="0.9" parameter="@clockRate.c:dplace" weight="0.1"/>
		ScaleOperator clockRateScaler = new ScaleOperator();
		clockRateScaler.initByName("scaleFactor", 0.9, "weight", 0.1, "parameter", clockRate, "upper", 1e10);

		List<Operator> operators = new ArrayList<>();
		operators.add(deltaExchanger);
		operators.add(clockRateScaler);
		//
		// <distribution id="posterior" spec="CompoundDistribution">
		// <distribution id="prior" spec="CompoundDistribution">
		// <distribution id="rates_prior" spec="beast.math.distributions.Prior"
		// x="@parameters">
		// <distr spec="tsa.parameterclone.helpers.RescaledDirichlet">
		// <!-- Use a flat distribution -->
		// <parameter name="sizes" idref="sizes" />
		// </distr>
		RescaledDirichlet rescaledDirichlet = new RescaledDirichlet(sizes);
		// </distribution>
		Prior ratesPrior = new Prior();
		ratesPrior.initByName("x", parameters, "distr", rescaledDirichlet);
		// </distribution>
		CompoundDistribution prior = new CompoundDistribution();
		prior.initByName("distribution", ratesPrior);
		prior.setID("prior");
		// <distribution id="likelihood" spec="CompoundDistribution">
		// <distribution id="treeSetLikelihood.dplace"
		// spec="beast.evolution.likelihood.ThreadedTreeSetLikelihood"
		// treeSetFile="/Users/remco/data/beast/global/scripts/gen/dplace-data/geo-rc/combined/geo-rc1197+almostnewwals+corrcals11.tree"
		// burnin="0">
		// <treeLikelihood id="treelk" spec="AncestralStateTreeLikelihood"
		// tree="@tree" data="@characters" tag="dplace">
		// <siteModel id="sitemodel" spec="SiteModel">
		// <substModel id="subst"
		// spec="tsa.correlatedcharacters.polycharacter.CorrelatedSubstitutionModel"
		// shape="2 2" alignment="@characters">
		// <rates id="actual_rates" spec="tsa.parameterclone.selector.Selector">
		// <groupings idref="indices" />
		// <parameters idref="parameters" />
		// <entry spec="beast.core.parameter.IntegerParameter">0 1 2 3 4 5 6 7
		// </entry>
		// </rates>
		Selector rates = new Selector(indices, parameters, new IntegerParameter("0 1 2 3 4 5 6 7"));
		// <frequencies name="frequencies" id="freqs"
		// spec="Frequencies" estimate="false" frequencies="@frequencies" />
		Frequencies freqs = new Frequencies();
		freqs.initByName("frequencies", frequencies);
		// </substModel>
		CorrelatedSubstitutionModel substModel = new CorrelatedSubstitutionModel(new IntegerParameter("2 2"),
				characters, rates, freqs);
		// </siteModel>
		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);
		// <branchRateModel id="StrictClockModel.c:dplace"
		// spec="beast.evolution.branchratemodel.StrictClockModel">
		// <parameter id="clockRate.c:dplace" estimate="true" name="clock.rate"
		// upper="1.0">0.5</parameter>
		// </branchRateModel>
		StrictClockModel clockModel = new StrictClockModel();
		clockModel.initByName("clock.rate", clockRate);
		// </treeLikelihood>
//		AncestralStateTreeLikelihood ancestralStateTreeLikelihood = new AncestralStateTreeLikelihood();
//		ancestralStateTreeLikelihood.initByName("branchRateModel", clockModel, "tree", tree, "data", characters, "tag",
//				"dplace", "siteModel", siteModel, "useJava", true);
		TreeLikelihood ancestralStateTreeLikelihood = new TreeLikelihood();
		ancestralStateTreeLikelihood.initByName("branchRateModel", clockModel, "tree", tree, "data", characters, "siteModel", siteModel);
		// </distribution>
		CompoundDistribution likelihood = new CompoundDistribution();
		likelihood.setID("likelihood");
		if (treeSetFileInput.get() != null) {
			ThreadedTreeSetLikelihood treeLikelihood = new ThreadedTreeSetLikelihood(ancestralStateTreeLikelihood, treeSetFileInput.get(), 0);
		// </distribution>
			likelihood.initByName("distribution", treeLikelihood);
		} else {
			likelihood.initByName("distribution", ancestralStateTreeLikelihood);
		}
		// </distribution>
		CompoundDistribution distribution = new CompoundDistribution();
		distribution.initByName("distribution", prior, "distribution", likelihood);
		distribution.setID("posterior");

		NS ns = new NS(5000 * particleCountInput.get(), 0, particleCountInput.get(), subChainLengthInput.get(), state, operators, distribution,
				epsilonInput.get());
		return ns;
	}

	private NS buildModelRelaxed(String indexValues) throws IOException {
		TreeParser tree = new TreeParser(null, newick, 0, false);

		Alignment data = processDataFile(tree.getTaxaNames());
		CompoundAlignment characters = new CompoundAlignment(data);

		RealParameter parameters = new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0");
		parameters.setBounds(0.05, 8.0);
		parameters.setID("parameters");

		IntegerParameter indices = new IntegerParameter(indexValues);
		indices.setID("indices");

		IntegerParameter sizes = new IntegerParameter("1 1 1 1 1 1 1 1");
		sizes.setID("sizes");

		RealParameter frequencies = new RealParameter("0.25 0.25 0.25 0.25");
		frequencies.setID("frequencies");

		RealParameter clockRate = new RealParameter("0.5");
		clockRate.setUpper(1.0);
		clockRate.setID("clockRate");

		RealParameter ucldStdev = new RealParameter("0.5");
		ucldStdev.setUpper(2.0);
		ucldStdev.setID("ucldStdev");
		
		IntegerParameter rateCategories = new IntegerParameter("1 1 1 1");
		rateCategories.setID("rateCategories");

		State state = new State();
		state.initByName("stateNode", parameters,
				"stateNode", ucldStdev,
				"stateNode", sizes,
				"stateNode",rateCategories,
				"stateNode", clockRate);

		DeltaExchangeOperator deltaExchanger = new DeltaExchangeOperator();
		deltaExchanger.initByName("weightvector", sizes, "parameter", parameters, "autoOptimize", false, "weight", 1.0);

		ScaleOperator clockRateScaler = new ScaleOperator();
		clockRateScaler.initByName("scaleFactor", 0.9, "weight", 0.1, "parameter", clockRate);


	    ScaleOperator stdDevScaler = new ScaleOperator();
		stdDevScaler.initByName("scaleFactor", 0.5, "weight", 0.1, "parameter", ucldStdev);		
		IntRandomWalkOperator categoriesRandomWalk = new IntRandomWalkOperator();
		categoriesRandomWalk.initByName("parameter", rateCategories, "weight", 1.0, "windowSize", 1);
		SwapOperator categoriesSwapOperator = new SwapOperator();
		categoriesSwapOperator.initByName("intparameter", rateCategories, "weight", 1.0);
		UniformOperator categoriesUniform = new UniformOperator();
		categoriesUniform.initByName("parameter", rateCategories, "weight", 1.0);
		
		List<Operator> operators = new ArrayList<>();
		operators.add(deltaExchanger);
		operators.add(clockRateScaler);
		operators.add(stdDevScaler);
		operators.add(categoriesRandomWalk);
		operators.add(categoriesSwapOperator);
		operators.add(categoriesUniform);


		RescaledDirichlet rescaledDirichlet = new RescaledDirichlet(sizes);

		Gamma gamma = new Gamma();
		gamma.initByName("alpha", "0.5396", "beta", "0.3819");
		Prior ucldStdevPrior = new Prior();
		ucldStdevPrior.initByName("x", ucldStdev, "distr", gamma);
		
		Prior ratesPrior = new Prior();		
		ratesPrior.initByName("x", parameters, "distr", rescaledDirichlet);

		CompoundDistribution prior = new CompoundDistribution();
		prior.initByName("distribution", ratesPrior, "distribution", ucldStdevPrior);
		prior.setID("prior");

		Selector rates = new Selector(indices, parameters, new IntegerParameter("0 1 2 3 4 5 6 7"));

		Frequencies freqs = new Frequencies();
		freqs.initByName("frequencies", frequencies);

		CorrelatedSubstitutionModel substModel = new CorrelatedSubstitutionModel(new IntegerParameter("2 2"),
				characters, rates, freqs);

		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);

		
		LogNormalDistributionModel lognormal = new LogNormalDistributionModel();
		lognormal.initByName("S", ucldStdev, "M", "1.0", "meanInRealSpace", true);		
		UCRelaxedClockModel clockModel = new UCRelaxedClockModel();
		clockModel.initByName("clock.rate", clockRate, "rateCategories", rateCategories, "tree", tree, "distr", lognormal);

		AncestralStateTreeLikelihood ancestralStateTreeLikelihood = new AncestralStateTreeLikelihood();
		ancestralStateTreeLikelihood.initByName("branchRateModel", clockModel, "tree", tree, "data", characters, "tag",
				"dplace", "siteModel", siteModel ,"useJava", true);

		CompoundDistribution likelihood = new CompoundDistribution();
		likelihood.setID("likelihood");
		if (treeSetFileInput.get() != null) {
			ThreadedTreeSetLikelihood treeLikelihood = new ThreadedTreeSetLikelihood(ancestralStateTreeLikelihood, treeSetFileInput.get(), 0);
			likelihood.initByName("distribution", treeLikelihood);
		} else {
			likelihood.initByName("distribution", ancestralStateTreeLikelihood);
		}

		CompoundDistribution distribution = new CompoundDistribution();
		distribution.initByName("distribution", prior, "distribution", likelihood);
		distribution.setID("posterior");

		NS ns = new NS(5000 * particleCountInput.get(), 0, particleCountInput.get(), subChainLengthInput.get(), state, operators, distribution,
				epsilonInput.get());
		return ns;
	}

	/**
	 * convert file exported from DPLACE into an Alignment
	 * 
	 * @param taxaNames
	 **/
	private Alignment processDataFile(String[] taxaNames) throws IOException {
		BufferedReader fin = new BufferedReader(new FileReader(dataFileInput.get()));
		String str = null;
		// skip first line
		str = fin.readLine();

		// determine columns containing the data
		str = fin.readLine();
		String[] strs = str.split(",");
		int socIndex = -1;
		int v1Input = -1;
		int v2Input = -1;
		String code = codeInput.get();
		for (int i = 0; i < strs.length; i++) {
			if (strs[i].equals(code)) {
				if (socIndex < 0) {
					socIndex = i;
				} else {
					fin.close();
					throw new IllegalArgumentException("more than one " + code + " column found");
				}
			}
			if (strs[i].startsWith("Code: ")) {
				if (v1Input < 0) {
					v1Input = i;
				} else if (v2Input < 0) {
					v2Input = i;
				} else {
					fin.close();
					throw new IllegalArgumentException("Found more than 2 Code: columns, but expected 2");
				}
			}

		}
		if (socIndex < 0) {
			fin.close();
			throw new IllegalArgumentException("Could not find " + codeInput + " column");
		}
		if (v1Input < 0) {
			fin.close();
			throw new IllegalArgumentException("Could not find any Code: column");
		}
		if (v2Input < 0) {
			fin.close();
			throw new IllegalArgumentException("Could only find one Code: column, but expected 2");
		}

		// get data out of columns
		Map<String, String> map = new LinkedHashMap<>();
		for (String iso : taxaNames) {
			map.put(iso, "??");
		}
		
		boolean [] values1 = toBoolean(thresholdV1Input.get());
		boolean [] values2 = toBoolean(thresholdV2Input.get());
		
		while (fin.ready()) {
			str = fin.readLine();
			// remove qouted text
			str = str.replaceAll("\"[^\"]+\"", "");
			strs = str.split(",");
			String iso = strs[socIndex];
			if (iso.length() > 0) {
				String v1org = strs[v1Input];
				String v1 = normalise(v1org, values1);
				String v2org = strs[v2Input];
				String v2 = normalise(v2org, values2);
				String valueOrg = v1org + v2org;
				String value = v1 + v2;
				// make sure that ambiguities are ignored
				if (v1.charAt(0)=='?' || v2.charAt(0)=='?') {
					value = "??";
				}
 				if (map.containsKey(iso)) {
					String value2 = map.get(iso);
					if (!value.equals(value2) && !value2.equals("??")) {
						Log.warning("Ambiguous value for " + iso + " (" + value + " & " + value2
								+ "). Setting value to \"??\"");
						map.put(iso, "??");
					} else {
						map.put(iso, value + " " + valueOrg);
					}
				}
			}

		}
		fin.close();

		// convert to Alignment
		Alignment data = new Alignment();
		data.userDataTypeInput.setValue(new Binary(), data);
		List<Sequence> seqs = data.sequenceInput.get();
		int [] siteCount = new int[4];
		for (String iso : map.keySet()) {
			String [] v = map.get(iso).split(" ");
			Sequence seq = new Sequence(iso, v[0]);
			if (v.length > 1) {
				seq.setID("seq" + seqs.size() + "_" + v[1]);
			} else {
				seq.setID("seq" + seqs.size() + "_??");  
			}
			seqs.add(seq);
			switch(map.get(iso)) {
			case "00" : siteCount[0]++;break;
			case "01" : siteCount[1]++;break;
			case "10" : siteCount[2]++;break;
			case "11" : siteCount[3]++;break;
			}
		}
		Log.info("value: count");
		Log.info("00   : " + siteCount[0]);
		Log.info("01   : " + siteCount[1]);
		Log.info("10   : " + siteCount[2]);
		Log.info("11   : " + siteCount[3]);
		this.siteCount = siteCount;
		
		data.initAndValidate();

		return data;
	}

	private boolean[] toBoolean(String str) {
		boolean [] b = new boolean[100];
		str = str.replaceAll("\\s", "");
		String [] strs = str.split(",");
		for (String s : strs) {
			b[Integer.parseInt(s)] = true;
		}
		return b;
	}

	private String normalise(String v, boolean [] threshold) {
		if (v.equals("NA")) {
			return "?";
		}
		try {
			int value = Integer.parseInt(v);
			if (threshold[value]) {
				return "1";
			} else {
				return "0";
			}
		} catch (NumberFormatException e) {
			return "?";
		}
	}

	public static void main(String[] args) throws Exception {
		new Application(new TSAModelSelector(), "DPLACE Model Selector", args);
		System.exit(0);
	}

}
