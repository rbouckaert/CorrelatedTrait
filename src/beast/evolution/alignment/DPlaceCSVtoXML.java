package beast.evolution.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.app.util.Application;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.Input.Validate;

public class DPlaceCSVtoXML extends Runnable {
	final public Input<File> fileInput = new Input<>("file", "file exported from DPlace database containing two features", Validate.REQUIRED);
	final public Input<String> code1Input = new Input<>("code1", "comma separated list of codes to be mapped to 1 for first feature", "1");
	final public Input<String> code2Input = new Input<>("code2", "comma separated list of codes to be mapped to 1 for second feature", "1");
	
	Map<String, String> data = new HashMap<>();

	static String [] isos = new String[]{"aaq","ach","aeb","agt","aka","ako","ale","aly","anq","aoc","apc","apd","apf","arq","ars","ary","axk","axl","bdy","bej","bem","bfa","biy","bkc","bla","bvr","cce","cea","ciw","col","coo","csi","ctc","dbl","ddj","dgc","dih","djb","dnj","drl","efe","esi","esk","ess","esu","eus","ewe","fij","fla","fln","fuv","gbd","gbu","gce","gia","gil","git","gup","guq","gur","gvn","gwj","gyd","haa","hax","hei","hgm","hoi","hrc","hur","ike","iow","jnj","kal","kat","kcg","kdd","kdj","kee","kej","khm","kjq","kju","kld","koy","ktw","lbz","len","lol","lrg","mah","mbr","mdj","mjd","mjp","mpb","mpj","mra","msh","mwp","naq","nay","ngc","ngh","nju","nmn","nol","nso","nsq","ntj","nuy","nya","nyy","ojb","ojs","oka","oma","ood","oon","orh","oun","pao","pbu","pcf","pes","piu","pjt","pko","plt","plw","por","ppi","ptw","ron","rxw","sbl","sge","shh","shi","shu","sii","skd","som","str","swh","tap","taq","tat","tew","thv","tig","tik","tjw","ttm","ttq","twf","tyu","udm","ukr","ump","unn","ute","vmv","wbt","wgu","wnw","woe","wrg","wro","wub","xeg","xet","xgf","xho","xnz","xtz","xuj","xuu","yak","yok","yor","yuf","yuq","zku","aar","abk","abr","abs","abt","abu","acm","acs","acv","ada","adi","adz","aec","aer","aes","afo","ahg","ahk","ahs","aim","ain","aiw","aji","ajp","akl","akp","alc","ald","aln","alz","amc","amf","amh","ami","anc","anu","any","aoi","aol","aon","apb","apj","apk","apl","apm","apn","apq","apw","ari","arn","arp","arv","arw","arz","asa","asb","ati","ats","atw","avi","awe","axb","ayl","ayn","ayr","bab","bak","bam","ban","bax","bbc","bbo","bbp","bca","bci","bcp","bcy","bdm","bdp","bdu","bea","bel","ben","bet","bex","bez","bfd","bfm","bfo","bfz","bgo","bgp","bhb","bho","bib","bif","bin","bip","biv","bja","bjg","bji","bjw","bkh","bkm","bko","bkq","bky","blc","ble","bmi","bmw","bnm","bnn","bnp","bod","bol","bom","bor","bot","boz","bpg","bpr","bqi","brf","bri","brm","brt","bsc","bsk","bst","bta","bud","buf","buk","bul","bun","buu","bux","bvb","bwd","bwr","bwu","bxk","bxr","byi","byn","byz","bza","bzd","bze","bzw","caa","cab","cad","caf","cak","cal","cao","caq","car","cax","cbi","ccg","cch","ccp","cde","cdr","ceb","ceg","ces","cgg","cha","chb","chc","chd","che","chh","chk","chl","cho","chp","chr","chv","chy","cic","cid","cjk","cjm","ckb","ckt","clc","clm","cmn","cni","coc","cod","com","cou","crd","crj","crk","cro","crq","crx","crz","csa","csm","csw","csz","ctm","cub","cup","cwe","dag","dak","daq","dav","dbw","dga","dgr","dhv","did","dif","dig","dik","dil","dim","diz","dje","dnt","dob","dos","dow","doz","dri","drs","dsh","dta","dtp","dua","dug","duj","dya","dyo","dyu","dzg","ebo","efi","ekg","ell","emp","ems","eng","enq","erg","est","ets","etu","evn","eya","eyo","fan","fip","fll","fmp","foi","fon","for","fra","fub","fuc","fud","fuf","fuh","fuj","fuy","fvr","gaa","gag","gax","gaz","gby","gde","ggu","giz","gld","gle","gnc","gno","goa","gog","gol","gow","grh","grt","gta","guc","guh","gui","guj","gun","guu","gux","guz","gwi","gya","gym","haq","has","hau","haw","hay","haz","hch","hdn","hea","heb","heh","her","hid","hke","hmb","hne","hnn","hoc","hop","hoz","hts","hun","hup","hus","huu","huv","hwo","hye","ian","iar","iba","ibb","ibo","ich","ida","idu","ifb","ife","igb","igl","iii","ijn","ikt","ilb","ing","irk","isl","iso","itl","its","izh","jac","jav","jbu","jib","jic","jiv","jle","jpn","kab","kac","kad","kai","kam","kao","kas","kay","kaz","kbc","kbd","kbk","kbl","kbp","kbr","kby","kca","kde","kdh","kdn","kdr","ked","kei","kel","ken","keo","ket","kfa","kff","kgo","kgq","kha","khk","khu","khy","kib","kic","kig","kij","kik","kin","kio","kjd","kki","kkz","kla","klb","klu","kmo","kmw","kna","knb","knk","kns","kod","kog","koo","koq","kor","kos","kot","kpc","kpg","kpj","kpy","kpz","kqi","kqn","kqq","kqs","kqy","krh","krn","kru","krw","ksb","ksf","ksv","kts","ktz","kui","kum","kun","kus","kut","kvb","kvj","kvn","kwd","kwe","kwk","kwy","kxc","kxu","kyh","kzy","lac","lag","laj","lak","lal","lam","lav","lbk","lbn","lch","lcp","lea","leb","led","leg","lel","lep","les","lgg","lgu","lic","lil","lin","lit","liv","lkr","lkt","lle","llx","lmi","lmp","lmu","lmw","lns","lnu","lob","log","lot","loz","lse","lua","lub","lue","lug","lui","lum","lun","luo","lus","maf","mai","mak","mal","mam","mas","mau","mav","maw","mbb","mbc","mch","mcn","mcu","mdd","mdi","mdp","mdy","mei","mek","men","mer","meu","mey","mez","mfq","mfz","mgd","mgo","mgp","mgr","mgs","mgu","mgz","mhi","mhj","mhk","mhq","mhr","mia","mic","mik","min","miq","mir","mjg","mji","mjw","mlq","mmg","mmk","mnc","mng","mnr","mns","mnv","moe","mor","mos","mot","mov","mox","mqu","mrc","mrh","mri","mrl","mrq","mrt","mrv","mrz","mtq","mtt","mua","mug","mus","mva","mvb","mvf","mvi","mvy","mwf","mwm","mwn","mwt","mwv","mxa","mya","mye","myh","myk","mym","myo","myu","myv","myx","mzb","mzm","mzv","nab","nak","nan","nap","nar","nau","nav","nbm","ncf","nci","ncm","ncz","ndc","nde","ndi","ndo","ndr","neb","neg","nez","ngb","ngi","ngo","ngp","nhn","nhr","nio","niq","niu","niv","njh","njm","njo","njy","nku","nld","nmg","nmu","nnw","nrb","nre","nsk","nsm","nsz","ntp","ntu","nuk","nup","nus","nwi","nyf","nyi","nyk","nym","nyn","nyo","nza","oaa","oac","ojc","ojg","ojv","ojw","oki","okv","old","olo","ona","opt","ori","oro","oss","ote","otr","otw","pab","pan","par","pau","paw","pbb","pbg","pbh","pbi","pbo","peb","pei","pej","pem","peq","pib","pid","pio","piw","pkb","pkp","pll","plu","pmt","pna","pnh","poc","pof","poi","pon","pot","pub","pwi","pwn","pwo","pyu","quc","qui","qun","quz","qxq","rad","ram","rap","rar","reg","res","rif","rim","rkh","rki","rtm","ruf","run","rus","rwm","rys","ryu","sac","sad","sah","saq","sat","sax","sbb","sbd","sbk","sbp","scl","scs","sda","sdj","sea","see","sef","seh","sei","sek","sel","ser","sez","sgc","sgw","shb","shc","she","shk","shp","shr","shs","sht","shy","sid","sil","sin","sis","siw","siz","sjw","skg","skt","slh","sln","slp","slu","smd","smj","smo","smw","sna","snd","sng","snk","snp","sns","sod","soo","sop","sot","soz","spa","squ","srp","srq","srr","srs","ssb","ssw","sub","suk","suq","sus","suw","sva","swp","syb","szb","taf","tah","tam","tao","tar","tau","tay","tbi","tbw","tbz","tca","tcc","tcx","tcz","tdx","teh","tel","tem","teo","ter","tet","tey","tfn","tha","thp","tht","tia","tic","til","tir","tiv","tiw","tix","tjm","tkg","tkl","tkm","tkp","tlb","tli","tll","tmj","tnq","tob","tod","tog","toi","tol","ton","toq","tos","tow","tpn","tpy","tpz","tqb","tqn","tsb","tsi","tsn","tso","tsv","tsz","ttj","ttr","ttv","tub","tuf","tui","tuk","tum","tuo","tuq","tur","tuu","tuv","tva","tvl","tvu","twa","twu","tzh","tzm","ubu","ude","ulc","uli","ulk","uma","umb","umo","unr","ure","uzn","vag","vai","ved","ven","vep","ver","vie","vka","vmw","vor","vot","vut","wac","wao","wap","was","waw","way","wba","wbp","wca","wic","wim","win","wiv","wiy","wlk","wln","wls","wlv","wnc","woc","wol","wrp","wuu","wuv","wya","wyb","xac","xal","xam","xav","xaw","xbr","xcw","xer","xnn","xog","xok","xom","xon","xpe","xsi","xsl","xsm","xsr","xsu","xvi","yad","yae","yaf","yag","yal","yao","yap","yaq","yar","yaz","yer","ykg","yle","ynn","yns","yom","yrk","yua","yuc","yue","yuk","yum","yup","yur","yyr","zen","zin","ziw","zkk","zlm","zmi","zmp","zms","zmz","zne","zpu","zul","zun"};
	
	@Override
	public void run() {
		List<String> tabu = new ArrayList<>();
		for (String str : new String [] {"san", "tgu", "wit", "djk", "daf", "afr", "dep", "srm", "egy", "kzh", "tmr", "hat", "lat"}) {
			tabu.add(str);
		}

		Set<String> code1 = new HashSet<>();
		String str = code1Input.get();
		String [] strs = str.split(",");
		for (String str2 : strs) {
			code1.add(str2.trim());
		}
		Set<String> code2 = new HashSet<>();
		str = code2Input.get();
		strs = str.split(",");
		for (String str2 : strs) {
			code2.add(str2.trim());
		}
		
		
		try {
			File file = fileInput.get();
	        BufferedReader fin = new BufferedReader(new FileReader(file));
	        str = null;
	        // eat up header
	        fin.readLine();
	        str = fin.readLine();
	        strs = AlignmentFromTraitForDPLace.split(str);
	        int col1 = -1;
	        int col2 = -1;
	        int iso = -1;
	        for (int i = 0; i < strs.length; i++) {
	        	if (strs[i].startsWith("Code")) {
	        		if (col1 >= 0) 
	        			col2 = i;
	        		else 
	        			col1 = i;
	        	}
	        	if (strs[i].matches("ISO.code")) {
	        		iso = i;
	        	}
	        }
	        if (col1 < 0 || col2 < 0 || iso < 0) {
	        	fin.close();
	        	throw new IllegalArgumentException("Could not find columns with Code: and/or ISO code in them. The file does not appear to be exported from DPLace");
	        }
	        while (fin.ready()) {
	            str = fin.readLine();
	            strs = AlignmentFromTraitForDPLace.split(str);
	            String s1 = strs[col1];
	            String s2 = strs[col2];
	            String iso0 = strs[iso];
	            if (iso0.length() >= 3) {
	            	if (!s1.equals("NA") && !s2.equals("NA")) {
		            	if (code1.contains(s1)) {
		            		if (code2.contains(s2)) {
		    	            	data.put(iso0, "11");
		            		} else {
		    	            	data.put(iso0, "10");
		            		}
		            	} else {
		            		if (code2.contains(s2)) {
		    	            	data.put(iso0, "01");
		            		} else {
		    	            	data.put(iso0, "00");
		            		}
		            	}
	            	} else {
    	            	data.put(iso0, "??");
	            		//tabu.add(iso0);
	            	}
	            }
	        }
	        fin.close();
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e.getMessage());
		}

		for (String iso0 : tabu) {
        	data.put(iso0, "??");
		}	
		
		for (String iso : isos) {
			String value = data.get(iso);
			if (value == null) {
				value = "??";
			}
			System.out.println("<sequence taxon=\"" + iso + "\" value=\"" + value+ "\"/>");
		}
	}

	@Override
	public void initAndValidate() {
	}
		
	public static void main(String[] args) throws Exception {
		new Application(new DPlaceCSVtoXML(), "Alignment From Trait or DPLace", args);
	}
}
