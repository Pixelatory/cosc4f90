package util;

import java.util.ArrayList;

public class Sequences {

    public static String[] basic1 = {
            "FFABCD",
            "ABCDFF",
            "GGABCD",
            "ABCDGG"
    };

    public static String[] basic2 = {
            "ABCDEF",
            "AXBXXCXXXDXEXXF",
            "ABCDEF"
    };

    public static String[] med2 = {
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
            "MAMBCMDEFGMHIJKLMNOMPQRMSTMUVWXYZ",
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    };

    public static String[] med3 = {
            "CBCADCAACE",
            "EACABDCADB",
            "DABAECBDCD",
            "DBEACEACCD",
            "DDABDEEEDE",
            "EEAECCAAEB",
            "EABEBCBCCB",
            "BAADDACDBB"
    };

    public static String[] spaces = {
            "AXXXXXBXXXXXCXXXXXDXXXXXEXXXXXFXXXXXGXXXXXHXXXXXI",
            "ZABCDERGHI",
            "PABCDERGHI"
    };

    /*
        Below here are sequences from balibase.
        Notes:
        - RV11 (very divergent, <20% identity)
        - RV12 (medium to divergent, 20-40% identity)
     */

    // In RV11/BB11006 (Actually called 1bbt_ac but cannot start var with number)
    public static String[] _1bbt_ac = {
            "GIFPVACSDGYGGLVTTDPKTADPVYGKVFNPPRNQLPGRFTNLLDVAEA",
            "CPTFLRFEGGVPYVTTKTDSDRVLAQFDMSLAAKHMSNTFLAGLAQYYTQ",
            "YSGTINLHFMFTGPTDAKARYMVAYAPPGMEPPKTPEAAAHCIHAEWDTG",
            "LNSKFTFSIPYLSAADYTYTASDVAETTNVQGWVCLFQITHGKADGDALV",
            "VLASAGKDFELRLPVDARAE"
    };

    // In RV12/BB12001
    public static String[] yua6_caeel = {
            "MVMSRTLVLVALLGFAYVCESALITNLPGAPISNFKQYSGYYNVGTKKNH",
            "MLHYWFVESQSNPSTDPVLLWLTGGPGCSGLSALLTEWGPWNVNTDGATL",
            "RTNPYSWNKNASILTLEAPAGVGYSYATDNNIATGDDQTASENWEALVAF",
            "FNEFPQYKGNDFYVTGESYGGIYVPTLVQTILDRQSQSHINIKGLAIGNG",
            "CVSANEGVDSLVNFLYHHGVVDQAKWEHMKTSCCHNDTDACPWHSFSEFS",
            "ACGEFVEATQQTAWNGGLNPYNMYADCISTSASFRFGMEYERRFNKKYTP",
            "EVLGTVPCLDESPVTNYLNRQDVRKALGIPSSLPAWSICSNAISYGYKRQ",
            "YGDMTSRVLNAVNNNNLKMMLYNGDVDLACNALMGQRFTDKLGLTLSKKK",
            "THFTVKGQIGGYVTQYKGSQVTFATVRGAGHMVPTDKPAVAEHIIQSFLF",
            "NKAF"
    };

    // In RV11/BB11002
    public static String[] labo_A = {
            "NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS",
            "NYITPVNS"
    };

    // In RV12/BB12009
    public static String[] CSPF_ECOLI = {
            "MSRKMTGIVKTFDGKSGKGLITPSDGRIDVQLHVSALNLRDAEEITTGLR",
            "VEFCRINGLRGPSAANVYLS"
    };

    // In RV12/BB12005
    public static String[] SODM_CANAL = {
            "MFSIRSSSRVLLKASSATTRATLNAAASKTFTRSKYSLPELDYEFSATEP",
            "YISGQINEIHYTKHHQTYVNNLNASIEQAVEAKSKGEVKKLVALEKAINF",
            "NGGGYLNHCLWWKNLAPVSQGGGQPPSEDSKLGKQIVKQFGSLDKLIEIT",
            "NGKLAGIQGSGWAFIVKNKANGDTIDVITTANQDTVTDPNLVPLIAIDAW",
            "EHAYYLQYQNVKADYFKNLWHVINWKEAERRFEF"
    };

    private static int numOfUniqueChars(String[] seq) {
        ArrayList<Character> uniqueChars = new ArrayList<>();
        int count = 0;

        for (String str : seq) {
            for (char c : str.toCharArray()) {
                if (!uniqueChars.contains(c)) {
                    uniqueChars.add(c);
                    count++;
                }
            }
        }

        return count;
    }

    public static void main(String[] args) {
        System.out.println(Sequences.numOfUniqueChars(Sequences._1bbt_ac));
        System.out.println(Sequences.numOfUniqueChars(Sequences.yua6_caeel));
        System.out.println(Sequences.numOfUniqueChars(Sequences.labo_A));
        System.out.println(Sequences.numOfUniqueChars(Sequences.CSPF_ECOLI));
        System.out.println(Sequences.numOfUniqueChars(Sequences.SODM_CANAL));
    }
}
