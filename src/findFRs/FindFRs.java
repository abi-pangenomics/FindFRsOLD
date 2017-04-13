package findFRs;

/**
 *
 * @author bmumey
 */
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import java.io.*;

public class FindFRs {

    static Graph g;
    static ArrayList<Sequence> sequences;
    static char[] fastaConcat;
    static int fastaConcatLen;
    static TreeMap<Integer, Integer> startToNode;
    static int[][] paths;

    static PriorityQueue<ClusterEdge> edgeQ;
    static PriorityQueue<ClusterNode> bestQ;
    static ClusterNode[] startNode;

    static ConcurrentHashMap<Integer, ClusterNode> nodeCluster;

    // command line arguments:
    static String dotFile = ""; // .dot filename
    static String fastaFile = ""; // .fasta filename
    static int K = -1; // k-mer size
    static double alpha = 0.0; // epsilon_r parameter
    //static int maxInsert; // maxInsert parameter
    static int minSup;
    static int minSize;
    //static String filePrefix = ""; // file name prefix

    static String[] colors = {"122,39,25", "92,227,60", "225,70,233", "100,198,222", "232,176,49", "50,39,85",
        "67,101,33", "222,142,186", "92,119,227", "206,225,151", "227,44,118", "229,66,41",
        "47,36,24", "225,167,130", "120,132,131", "104,232,178", "158,43,133", "228,228,42", "213,217,213",
        "118,64,79", "88,155,219", "226,118,222", "146,197,53", "222,100,89", "224,117,41", "160,96,228",
        "137,89,151", "126,209,119", "145,109,70", "91,176,164", "54,81,103", "164,174,137", "172,166,48",
        "56,86,143", "210,184,226", "175,123,35", "129,161,88", "158,47,85", "87,231,225", "216,189,112", "49,111,75",
        "89,137,168", "209,118,134", "33,63,44", "166,128,142", "53,137,55", "80,76,161", "170,124,221", "57,62,13",
        "176,40,40", "94,179,129", "71,176,51", "223,62,170", "78,25,30", "148,69,172", "122,105,31", "56,33,53",
        "112,150,40", "239,111,176", "96,55,25", "107,90,87", "164,74,28", "171,198,226", "152,131,176", "166,225,211",
        "53,121,117", "220,58,86", "86,18,56", "225,197,171", "139,142,217", "216,151,223", "97,229,117", "225,155,85",
        "31,48,58", "160,146,88", "185,71,129", "164,233,55", "234,171,187", "110,97,125", "177,169,175", "177,104,68",
        "97,48,122", "237,139,128", "187,96,166", "225,90,127", "97,92,55", "124,35,99", "210,64,194", "154,88,84",
        "100,63,100", "140,42,54", "105,132,99", "186,227,103", "224,222,81", "191,140,126", "200,230,182", "166,87,123",
        "72,74,58", "212,222,124", "205,52,136"};
    static String[] svgcolors = {"aliceblue", "antiquewhite", "aqua", "aquamarine",
        "azure", "beige", "bisque", "black", "blanchedalmond", "blue",
        "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse",
        "chocolate", "coral", "cornflowerblue", "cornsilk", "crimson",
        "cyan", "darkblue", "darkcyan", "darkgoldenrod", "darkgray",
        "darkgreen", "darkgrey", "darkkhaki", "darkmagenta", "darkolivegreen",
        "darkorange", "darkorchid", "darkred", "darksalmon", "darkseagreen",
        "darkslateblue", "darkslategray", "darkslategrey", "darkturquoise", "darkviolet",
        "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue",
        "firebrick", "floralwhite", "forestgreen", "fuchsia", "gainsboro",
        "ghostwhite", "gold", "goldenrod", "gray", "grey",
        "green", "greenyellow", "honeydew", "hotpink", "indianred",
        "indigo", "ivory", "khaki", "lavender", "lavenderblush",
        "lawngreen", "lemonchiffon", "lightblue", "lightcoral", "lightcyan",
        "lightgoldenrodyellow", "lightgray", "lightgreen", "lightgrey", "lightpink",
        "lightsalmon", "lightseagreen", "lightskyblue", "lightslategray", "lightslategrey",
        "lightsteelblue", "lightyellow", "lime", "limegreen", "linen",
        "magenta", "maroon", "mediumaquamarine", "mediumblue", "mediumorchid",
        "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise",
        "mediumvioletred", "midnightblue", "mintcream", "mistyrose", "moccasin",
        "navajowhite", "navy", "oldlace", "olive", "olivedrab",
        "orange", "orangered", "orchid", "palegoldenrod", "palegreen",
        "paleturquoise", "palevioletred", "papayawhip", "peachpuff", "peru",
        "pink", "plum", "powderblue", "purple", "red",
        "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown",
        "seagreen", "seashell", "sienna", "silver", "skyblue",
        "slateblue", "slategray", "slategrey", "snow", "springgreen",
        "steelblue", "tan", "teal", "thistle", "tomato",
        "turquoise", "violet", "wheat", "white", "whitesmoke",
        "yellow", "yellowgreen"};

    static void readData() {
        g = ReadInput.readDotFile(dotFile);

        startToNode = new TreeMap<Integer, Integer>();
        for (int i = 0; i < g.starts.length; i++) {
            for (int j = 0; j < g.starts[i].length; j++) {
                startToNode.put(g.starts[i][j], i);
            }
            int firstStart = g.starts[i][0];
            g.starts[i] = new int[1]; // only save 1st start
            g.starts[i][0] = firstStart;
        }

        sequences = ReadInput.readFastaFile(fastaFile);
    }

//    static int[] removeDuplicates(int[] a) {
//        int[] b = null;
//
//        if (a.length > 0) {
//            ArrayList<Integer> B = new ArrayList<Integer>();
//            B.add(a[0]);
//            for (int i = 1; i < a.length; i++) {
//                if (a[i] != B.get(B.size() - 1)) {
//                    B.add(a[i]);
//                }
//            }
//            int i = 0;
//            b = new int[B.size()];
//            for (Integer pobj : B) {
//                b[i++] = pobj;
//            }
//        }
//        return b;
//    }
    static void buildPaths() {
        ArrayList<ArrayList<Integer>> pathsAL = new ArrayList<ArrayList<Integer>>();
        int curStart = 1;
        int seqStart = 1;
        int seqEnd;
        int prevStart = 0;
        for (Sequence s : sequences) {
            s.startPos = seqStart;
            s.length = s.seq.length();
            seqEnd = seqStart + s.length - 1;
            curStart = seqStart;
            boolean seqendCovered = false;
            ArrayList path = new ArrayList<Integer>();
            while (!seqendCovered) {
                path.add(startToNode.get(curStart));
                if (curStart + g.length[startToNode.get(curStart)] - 1 >= seqEnd) {
                    seqendCovered = true;
                }
                curStart += g.length[startToNode.get(curStart)] - (K - 1);

            }
            pathsAL.add(path);
            //System.out.println(path);
            seqStart = seqEnd + 2;

            fastaConcatLen += 1 + s.seq.length();
        }
        fastaConcatLen++;
        fastaConcat = new char[fastaConcatLen];
        for (Sequence s : sequences) {
            fastaConcat[s.startPos - 1] = '$';
            for (int i = 0; i < s.length; i++) {
                fastaConcat[s.startPos + i] = s.seq.charAt(i);
            }
            s.seq = null; // no longer needed
        }
        fastaConcat[fastaConcat.length - 1] = '$';
        //System.out.println(Arrays.toString(fastaConcat));
        System.out.println("number of paths: " + pathsAL.size());

        paths = new int[pathsAL.size()][];
        for (int i = 0; i < pathsAL.size(); i++) {
            ArrayList<Integer> path = pathsAL.get(i);
            paths[i] = new int[path.size()];
            for (int j = 0; j < paths[i].length; j++) {
                paths[i][j] = path.get(j);
            }
        }

        pathsAL.clear();
        pathsAL = null; // can be gc'ed

        System.out.println("finding node paths");

        g.containsN = new boolean[g.numNodes];
        for (int i = 0; i < g.numNodes; i++) {
            g.containsN[i] = false;
            for (int j = 0; j < g.length[i]; j++) {
                if (fastaConcat[g.starts[i][0] + j] == 'N') {
                    g.containsN[i] = true;
                    break;
                }
            }
        }

        // find paths for each node:
        g.nodePaths = new TreeMap<Integer, TreeSet<Integer>>();
        for (int i = 0; i < g.numNodes; i++) {
            if (!g.containsN[i]) {
                g.nodePaths.put(i, new TreeSet<Integer>());
            }
        }

        for (int i = 0; i < paths.length; i++) {
            for (int j = 0; j < paths[i].length; j++) {
                if (!g.containsN[paths[i][j]]) {
                    g.nodePaths.get(paths[i][j]).add(i);
                }
            }
        }

    }

    static int gap(int[] path, int start, int stop) {
        int curStartLoc = 1;
        int curEndLoc = 1;
        for (int i = start; i <= stop; i++) {
            curEndLoc = curStartLoc + g.length[path[i]] - 1;
            curStartLoc += g.length[path[i]] - (K - 1);
        }
        return Math.max(0, curEndLoc - g.length[path[start]] - g.length[path[stop]]);
    }

    static ConcurrentLinkedQueue<PathSegment> computeSupport(ClusterNode clust, boolean createPSList) {
        if (clust.pathLocs == null) {
            clust.findPathLocs();
        }
        ConcurrentLinkedQueue<PathSegment> segList = new ConcurrentLinkedQueue<PathSegment>();
        AtomicInteger clustSupport = new AtomicInteger(0);
        clust.pathLocs.keySet().parallelStream().forEach((P) -> {
            int[] locs = clust.pathLocs.get(P);
            int start = 0;
            while (start < locs.length) {
                int last = start;
                while (last + 1 < locs.length && locs[last + 1] == locs[last] + 1) {
                    last++;
                }
                if (locs[last] - locs[start] + 1 >= alpha * clust.size) {
                    clustSupport.incrementAndGet();
                    if (createPSList) {
                        segList.add(new PathSegment(P, locs[start], locs[last]));
                    }
                }
                start = last + 1;
            }
        });
        clust.support = clustSupport.get();
        if (createPSList) {
            return segList;
        }
        return null;
    }

//    static ConcurrentLinkedQueue<PathSegment> findSupport(ClusterNode clust) {
//        ConcurrentLinkedQueue<PathSegment> segList = new ConcurrentLinkedQueue<PathSegment>();
//        ConcurrentHashMap<Integer, TreeSet<Integer>> cPL;
//        if (clust.pathLocs != null) {
//            cPL = clust.pathLocs;
//        } else {
//            cPL = new ConcurrentHashMap<Integer, TreeSet<Integer>>();
//            clust.addLocs(cPL);
//            clust.pathLocs = cPL;
//        }
//        cPL.keySet().parallelStream().forEach((P) -> {
//            TreeSet<Integer> locs = cPL.get(P);
//            Iterator<Integer> iter = locs.iterator();
//            Integer start = iter.next();
//            Integer last = start;
//            while (iter.hasNext()) {
//                Integer next = iter.next();
//                if (next > last + 1) {
//                    if (last - start + 1 >= alpha * clust.size) {
//                        PathSegment ps = new PathSegment();
//                        ps.path = P;
//                        ps.start = start;
//                        ps.stop = last;
//                        segList.add(ps);
//                    }
//                    start = next;
//                }
//                last = next;
//            }
//            if (last - start + 1 >= alpha * clust.size) {
//                PathSegment ps = new PathSegment();
//                ps.path = P;
//                ps.start = start;
//                ps.stop = last;
//                segList.add(ps);
//            }
//        });
//        return segList;
//    }
//    static void addToQueue(ClusterNode clust) {
//        clust.size = clust.left.size + clust.right.size;
//        computeSupport(clust, false);
//        if (clust.support >= 0.9 * minSup) {
//            edgeQ.add(clust);
//        }
//    }
    static void finalizeEdge(ClusterEdge e) {
        ClusterNode newRoot = new ClusterNode();
        newRoot.left = e.u;
        newRoot.right = e.v;
        newRoot.parent = null;
        newRoot.pathLocs = null;
        newRoot.edges = new ArrayList<ClusterEdge>();
        newRoot.size = newRoot.left.size + newRoot.right.size;

        newRoot.left.parent = newRoot;
        newRoot.right.parent = newRoot;
        computeSupport(newRoot, false);
//        if (newRoot.left.pathLocs != null && newRoot.left.node == -1) {
//            newRoot.left.pathLocs.clear();
//            newRoot.left.pathLocs = null;
//        }
//        if (newRoot.right.pathLocs != null && newRoot.right.node == -1) {
//            newRoot.right.pathLocs.clear();
//            newRoot.right.pathLocs = null;
//        }

        TreeSet<ClusterNode> neighbors = new TreeSet<ClusterNode>();
        for (int i = 0; i < e.u.edges.size(); i++) {
            ClusterEdge ce = e.u.edges.get(i);
            ClusterNode n = ce.other(e.u);
            if (n != e.v && n.parent == null) {
                neighbors.add(n);
            }
            edgeQ.remove(ce);
            n.edges.remove(ce);
        }
        for (int i = 0; i < e.v.edges.size(); i++) {
            ClusterEdge ce = e.v.edges.get(i);
            ClusterNode n = ce.other(e.v);
            if (n != e.u && n.parent == null) {
                neighbors.add(n);
            }
            edgeQ.remove(ce);
            n.edges.remove(ce);
        }
        for (ClusterNode n : neighbors) {
            ClusterNode tmpClst = new ClusterNode();
            tmpClst.left = newRoot;
            tmpClst.right = n;
            tmpClst.parent = null;
            tmpClst.size = tmpClst.left.size + tmpClst.right.size;
            computeSupport(tmpClst, false);
            tmpClst.pathLocs.clear();
            ClusterEdge newE = new ClusterEdge(tmpClst.left, tmpClst.right, tmpClst.support);
            newRoot.edges.add(newE);
            n.edges.add(newE);
            edgeQ.add(newE);
        }
    }

    static void exploreSolns() {
        System.out.println("creating node clusters");
        nodeCluster = new ConcurrentHashMap<Integer, ClusterNode>(g.numNodes);

        // create initial node clusters
        g.nodePaths.keySet().parallelStream().forEach((N) -> {
            ClusterNode nodeClst = new ClusterNode();
            nodeClst.parent = nodeClst.left = nodeClst.right = null;
            nodeClst.node = N;
            nodeClst.size = 1;
            nodeClst.pathLocs = new ConcurrentHashMap<Integer, int[]>();
            nodeClst.edges = new ArrayList<ClusterEdge>();
            nodeCluster.put(N, nodeClst);
        });
        g.nodePaths.clear(); //not used after this
        g.nodePaths = null;
        TreeMap<Integer, TreeMap<Integer, ArrayList<Integer>>> nodePathLocs
                = new TreeMap<Integer, TreeMap<Integer, ArrayList<Integer>>>();
        for (int p = 0; p < paths.length; p++) {
            for (int i = 0; i < paths[p].length; i++) {
                int n = paths[p][i];
                if (!nodePathLocs.containsKey(n)) {
                    nodePathLocs.put(n, new TreeMap<Integer, ArrayList<Integer>>());
                }
                if (!nodePathLocs.get(n).containsKey(p)) {
                    nodePathLocs.get(n).put(p, new ArrayList<Integer>());
                }
                nodePathLocs.get(n).get(p).add(i);
            }
        }
        for (Integer node : nodePathLocs.keySet()) {
            nodePathLocs.get(node).keySet().parallelStream().forEach((P) -> {
                ArrayList<Integer> al = nodePathLocs.get(node).get(P);
                int[] ar = new int[al.size()];
                int i = 0;
                for (Integer x : al) {
                    ar[i++] = x;
                }
                if (nodeCluster.containsKey(node)) {
                    nodeCluster.get(node).pathLocs.put(P, ar);
                }
            });
        }
        nodePathLocs.clear();

        System.out.println("computing node support");
        for (ClusterNode c : nodeCluster.values()) {
            computeSupport(c, false);
        }
        // create initial edges
        edgeQ = new PriorityQueue<ClusterEdge>();
        System.out.println("adding edges");
        for (Integer N : nodeCluster.keySet()) {
            for (int i = 0; i < g.neighbor[N].length; i++) {
                if (nodeCluster.containsKey(g.neighbor[N][i])) {
                    ClusterNode tmpClst = new ClusterNode();
                    tmpClst.left = nodeCluster.get(N);
                    tmpClst.right = nodeCluster.get(g.neighbor[N][i]);
                    tmpClst.parent = null;
                    tmpClst.size = 2;
                    computeSupport(tmpClst, false);
                    tmpClst.pathLocs.clear();
                    ClusterEdge newE = new ClusterEdge(tmpClst.left, tmpClst.right, tmpClst.support);
                    edgeQ.add(newE);
                }
            }
        }
        System.out.println("edge queue size: " + edgeQ.size());
        ClusterEdge e;
        int count = 0;
        while ((e = edgeQ.poll()) != null) {
            if (e.potentialSup > 0 && e.u.parent == null && e.v.parent == null) {
                finalizeEdge(e);
            }
            count++;
            if (count % 1000 == 0) {
                System.out.println("# finalized: " + count);
            }
        }
        System.out.println("finding root FRs");
        HashSet<ClusterNode> roots = new HashSet<ClusterNode>();
        for (ClusterNode leaf : nodeCluster.values()) {
            ClusterNode cur = leaf;
            while (cur.parent != null) {
                cur = cur.parent;
            }
            roots.add(cur);
        }
        bestQ = new PriorityQueue<ClusterNode>();
        System.out.println("number of root FRs: " + roots.size());
        for (ClusterNode root : roots) {
            findBestSolns(root, 0);
        }
        System.out.println("number of iFRs: " + bestQ.size());
        outputResults();
    }

    static void findBestSolns(ClusterNode clust, int parentSup) {
        if (clust.edges != null) {
            clust.edges.clear();
            clust.edges = null;
        }
        if (clust.support > parentSup && clust.size >= minSize && clust.support >= minSup) {
            bestQ.add(clust);
        }
        if (clust.left != null) {
            findBestSolns(clust.left, Math.max(parentSup, clust.support));
        }
        if (clust.right != null) {
            findBestSolns(clust.right, Math.max(parentSup, clust.support));
        }
    }

    static int[] findFastaLoc(PathSegment ps) {
        int[] startStop = new int[2];
        int curPos = sequences.get(ps.path).startPos;

        int curIndex = 0;
        while (curIndex != ps.start) {
            curPos += g.length[startToNode.get(curPos)] - (K - 1);
            curIndex++;
        }
        startStop[0] = curPos - sequences.get(ps.path).startPos; // assume fasta seq indices start at 0
        while (curIndex != ps.stop) {
            curPos += g.length[startToNode.get(curPos)] - (K - 1);
            curIndex++;
        }
        startStop[1] = curPos + g.length[startToNode.get(curPos)] - 1 - sequences.get(ps.path).startPos + 1; // last position is excluded in BED format
        return startStop; // *** TODO *** deal with $ in final kmer node
    }

//    static void printCluster(ClusterNode cluster) {
//        ConcurrentLinkedQueue<PathSegment> supportingSegments = findSupport(cluster);
//        System.out.println("cluster nodes: " + cluster.getNodeSet());
//        System.out.println("total support: " + cluster.support);
//
//        for (PathSegment ps : supportingSegments) {
//            System.out.println("from fasta seq: " + sequences.get(ps.path).label);
//            System.out.println("node subpath: ");
//            for (int i = ps.start; i <= ps.stop; i++) {
//                System.out.print(" " + paths[ps.path][i]);
//            }
//            System.out.println();
//            System.out.print("fasta location: ");
//            int[] startStop = findFastaLoc(ps);
//            System.out.println(startStop[0] + "," + startStop[1]);
//        }
//        System.out.println();
//    }
    int findLen(int[] pa) {
        int c = 0;
        for (int i = 0; i < pa.length; i++) {
            if (i == 0) {
                c += g.length[pa[i]];

            } else {
                c += g.length[pa[i]] - (K - 1);
            }
        }
        return c;
    }

//    static String getStrainForSeq(int i) {
//        String[] tmp = sequences.get(i).label.split("_");
//        String strain;
//        if (tmp[1].startsWith("gi")) {
//            strain = tmp[0]; // for 10 yeast
//        } else {
//            strain = tmp[0] + "_" + tmp[1] + "_" + tmp[2]; // for 55 yeast
//        }
//        return strain;
//    }
    static void clearCPL(ClusterNode clust) {
        if (clust.node == -1) {
            if (clust.pathLocs != null) {
                clust.pathLocs.clear();
            }
            clust.pathLocs = null;
        }
        if (clust.left != null) {
            clearCPL(clust.left);
        }
        if (clust.right != null) {
            clearCPL(clust.right);
        }
    }

    static void outputResults() {
        ClusterNode top;
        ArrayList<ClusterNode> iFRs = new ArrayList<ClusterNode>();
        while ((top = bestQ.poll()) != null) {
            iFRs.add(top);
        }

        try {
            String paramString = "-k" + K + "-a" + alpha + "-sup" + minSup + "-size" + minSize;
            String[] tmp = dotFile.split("/");
            String dotName = tmp[tmp.length - 1];
            tmp = fastaFile.split("/");
            String fastaName = tmp[tmp.length - 1];
            String filePrefix = dotName + "-" + fastaName;
            String rd = "results-" + filePrefix + "/";
            File resultsDir = new File(rd);
            resultsDir.mkdir();

            HashMap<Integer, TreeSet<Integer>> nodeFRset = new HashMap<Integer, TreeSet<Integer>>();
            BufferedWriter frOut = new BufferedWriter(new FileWriter(rd + filePrefix + paramString + ".frs.txt"));
            for (int fr = 0; fr < iFRs.size(); fr++) {
                ClusterNode iFR = iFRs.get(fr);
                String frName = "fr-" + fr;
                TreeSet<Integer> clustNodes = iFR.getNodeSet();
                frOut.write(frName);
                for (Integer n : clustNodes) {
                    frOut.write("," + n);
                    if (!nodeFRset.containsKey(n)) {
                        nodeFRset.put(n, new TreeSet<Integer>());
                    }
                    nodeFRset.get(n).add(fr);
                }
                frOut.write("\n");
            }
            frOut.close();

            System.out.println("writing bed file");
            BufferedWriter bedOut = new BufferedWriter(new FileWriter(rd + filePrefix + paramString + ".bed"));
            TreeMap<String, TreeMap<Integer, Integer>> seqFRcount = new TreeMap<String, TreeMap<Integer, Integer>>();
            TreeMap<String, TreeMap<Integer, LinkedList<String>>> seqIndxFRstr = new TreeMap<String, TreeMap<Integer, LinkedList<String>>>();
            TreeMap<Integer, Integer> frAvgLen = new TreeMap<Integer, Integer>();
            for (int fr = 0; fr < iFRs.size(); fr++) {
                ClusterNode iFR = iFRs.get(fr);
                if ((fr % 100) == 0) {
                    System.out.println("writing fr-" + fr);
                }
                ConcurrentLinkedQueue<PathSegment> supportingSegments = computeSupport(iFR, true);
                //System.out.println("finished computing support");
//                if (iFR.parent == null) {
//                    clearCPL(iFR);
//                }
                frAvgLen.put(fr, 0);
                for (PathSegment ps : supportingSegments) {
                    String name = sequences.get(ps.path).label;
                    if (!seqFRcount.containsKey(name)) {
                        seqFRcount.put(name, new TreeMap<Integer, Integer>());
                    }
                    if (!seqFRcount.get(name).containsKey(fr)) {
                        seqFRcount.get(name).put(fr, 0);
                    }
                    seqFRcount.get(name).put(fr, seqFRcount.get(name).get(fr) + 1);
                    int[] startStop = findFastaLoc(ps);
                    int frLen = startStop[1] - startStop[0]; // last position is excluded
                    frAvgLen.put(fr, frLen + frAvgLen.get(fr));

                    bedOut.write(name // chrom
                            + "\t" + startStop[0] // chromStart (starts with 0)
                            + "\t" + startStop[1] // chromEnd
                            + "\t" + "fr-" + fr// name
                            + "\t" + Math.round(iFR.support) // score
                            + "\t+" // strand
                            + "\t" + 0 // thickstart
                            + "\t" + 0 // thickend
                            + "\t" + colors[fr % colors.length] // itemRGB
                            + "\t" + frLen // FR length
                            + "\n");
                    if (!seqIndxFRstr.containsKey(name)) {
                        seqIndxFRstr.put(name, new TreeMap<Integer, LinkedList<String>>());
                    }
                    if (!seqIndxFRstr.get(name).containsKey(ps.start)) {
                        seqIndxFRstr.get(name).put(ps.start, new LinkedList<String>());
                    }
                    if (!seqIndxFRstr.get(name).containsKey(ps.stop)) {
                        seqIndxFRstr.get(name).put(ps.stop, new LinkedList<String>());
                    }
                    seqIndxFRstr.get(name).get(ps.start).addFirst(" [fr-" + fr + ":" + startStop[0]);
                    seqIndxFRstr.get(name).get(ps.stop).addLast(" fr-" + fr + ":" + startStop[1] + "] ");
                }
                frAvgLen.put(fr, frAvgLen.get(fr) / iFR.support);
            }
            bedOut.close();

            System.out.println("writing dist file");
            BufferedWriter distOut = new BufferedWriter(new FileWriter(rd + filePrefix + paramString + ".dist.txt"));
            distOut.write("FR,size,support,avg length\n");
            for (int fr = 0; fr < iFRs.size(); fr++) {
                ClusterNode iFR = iFRs.get(fr);
                distOut.write("fr-" + fr + "," + iFR.size + "," + iFR.support + "," + frAvgLen.get(fr) + "\n");
            }
            distOut.close();

            System.out.println("writing frpaths file");
            BufferedWriter frPathsOut = new BufferedWriter(new FileWriter(rd + filePrefix + paramString + ".frpaths.txt"));
            for (int i = 0; i < paths.length; i++) {
                String name = sequences.get(i).label;
                if (seqIndxFRstr.containsKey(name)) {
                    frPathsOut.write(name + ",");
                    for (int pos : seqIndxFRstr.get(name).keySet()) {
                        LinkedList<String> ll = seqIndxFRstr.get(name).get(pos);
                        for (String s : ll) {
                            frPathsOut.write(s);
                        }
                    }
                    frPathsOut.write("\n");
                }
            }
            frPathsOut.close();

            System.out.println("writing csfr file");
            BufferedWriter seqFROut = new BufferedWriter(new FileWriter(rd + filePrefix + paramString + ".csfr.txt"));
            for (String seq : seqFRcount.keySet()) {
                seqFROut.write(seq);
                for (Integer F : seqFRcount.get(seq).keySet()) {
                    seqFROut.write("," + F + ":" + seqFRcount.get(seq).get(F));
                }
                seqFROut.write("\n");
            }
            seqFROut.close();

            System.out.println("done");
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(-1);
        }

    }

    public static void main(String[] args) {
        // parse args:
        if (args.length != 6) {
            System.out.println("Usage: java findFRs dotFile faFile K alpha minsup minsize");
            System.out.println(Arrays.toString(args));
            System.exit(0);
        }
        dotFile = args[0];
        fastaFile = args[1];
        K = Integer.parseInt(args[2]);
        alpha = Double.parseDouble(args[3]);
        minSup = Integer.parseInt(args[4]);
        minSize = Integer.parseInt(args[5]);
        //maxInsert = Integer.parseInt(args[2]);

        readData();
        buildPaths();

        startNode = new ClusterNode[g.numNodes];
        exploreSolns();
    }
}
