package phylolab.mrp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import phylolab.NewickTokenizer;

public class FastMRP {

    private static final String NEXUS = "NEXUS";
    private static final String NEXUS_HEADER = "#NEXUS\n"
            + "begin data;\n"
            + "\t dimensions ntax = @ nchar = @;\n"
            + "\t format missing = ?;\n"
            + "\tmatrix\n";
    private static final String NEXUS_FOOTER = "\t;\n"
            + "end;\n";
    private static final String PHYLIP = "PHYLIP";
    private static final String PHYLIP_HEADER = "@ @\n";
    private final String TREES_FILENAME;
    private final String MRP_FILENAME;
    private final String FORMAT;
    private Random random;
    private int treeCount;
    private final InfoPerTaxa PER_TAXA_INFO;
    private final TreeEndIndex TREE_END_INDEX;
    /** First list is the STACK position, each Collection is a bipartition 
     * (i.e index of taxa in one of the partitions) */
    private final LinkedList<Collection<Integer>> STACK;
    private char one;
    private char zero;
    private char missing;

    public FastMRP(
            final String IN_FILENAME, 
            final String OUT_FILENAME, 
            final String FORMAT) {
        this.TREES_FILENAME = IN_FILENAME;
        this.MRP_FILENAME = OUT_FILENAME;
        this.FORMAT = FORMAT;
        this.PER_TAXA_INFO = new InfoPerTaxa();
        this.TREE_END_INDEX = new TreeEndIndex();
        this.STACK = new LinkedList<>();
        this.one = '1';
        this.zero = '0';
        this.missing = '?';
    }

    public void setCharacters(char newOne, char newZero, char newMissing) {
        one = newOne;
        zero = newZero;
        missing = newMissing;
    }

    private void readTreesFile() throws IOException {
        BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(TREES_FILENAME)));
        String tree;
        int treeId = 0;
        int lastColumnInd = -1;
        //Read File Line By Line
        while ((tree = in.readLine()) != null) {
            NewickTokenizer tokenizer = new NewickTokenizer(tree);
            if (!"(".equals(tokenizer.nextToken())) {
                throw new RuntimeException("The tree does not start with a (");
            }
            while (tokenizer.hasNext()) {
                String token = tokenizer.nextToken();
                
                if (token != null) switch (token) {
                    case "(":
                        STACK.addLast(new ArrayList<>());
                        break;
                    case ")":
                        if (STACK.size() > 0) {
                            Collection<Integer> top = STACK.getLast();
                            lastColumnInd = PER_TAXA_INFO.addBipartition(top);
                            STACK.removeLast();
                            if (STACK.size() > 0) {
                                STACK.getLast().addAll(top);
                            }
                        }
                        break;
                    case ";":
                        TREE_END_INDEX.addEndIndex(lastColumnInd);
                        break;
                    default:
                        Integer seqId = PER_TAXA_INFO.mapNamesToIds(token);
                        if (STACK.size() > 0) {
                            STACK.getLast().add(seqId);
                        }
                        PER_TAXA_INFO.addTreeToSequence(seqId, treeId);
                        break;
                }
            }
            treeId++;
        }
        treeCount = treeId;
    }

    private void writeMRPToFile() throws IOException {
        try (BufferedWriter OUT = new BufferedWriter(
                new FileWriter(MRP_FILENAME))) {
            writeHeaderToFile(OUT);
            // iterate over sequences, one by one
            PER_TAXA_INFO.startSeqIteration();
            while (PER_TAXA_INFO.hasMoreTaxa()) {
                Integer nextOneColumn;
                
                PER_TAXA_INFO.nextSequence();
                writeSequenceNameToFile(OUT);
                nextOneColumn = 
                        PER_TAXA_INFO.nextBipartitionIndexForCurrentSequence();
                TREE_END_INDEX.restartTreeIteration();
                for (int tree = 0, column = 0; tree < treeCount; tree++) {
                    int treeEndIndex = TREE_END_INDEX.getNexTreeEndInd();
                    
                    if (PER_TAXA_INFO.currentSeqIsInTree(tree)) {
                        while (column <= treeEndIndex) {
                            Boolean coding = 
                                    PER_TAXA_INFO.columnCoding.get(column);
                            
                            if (nextOneColumn == column) {
                                OUT.write(coding ? one : zero);
                                nextOneColumn = PER_TAXA_INFO.nextBipartitionIndexForCurrentSequence();
                            } else {
                                OUT.write(coding ? zero : one);
                            }
                            column++;
                        }
                    } else {
                        char[] missings = new char[treeEndIndex - column + 1];
                        
                        Arrays.fill(missings, missing);
                        OUT.write(missings);
                        column = treeEndIndex + 1;
                    }
                }
                OUT.newLine();
            }
            writeFooterToFile(OUT);
            OUT.flush();
        }
    }

    private void writeFooterToFile(BufferedWriter out) throws IOException {
        if (NEXUS.equalsIgnoreCase(FORMAT)) {
            out.write(NEXUS_FOOTER);
        }
    }

    private void writeHeaderToFile(BufferedWriter out) throws IOException {
        if (NEXUS.equalsIgnoreCase(FORMAT)) {
            String header = NEXUS_HEADER.replaceFirst("@", PER_TAXA_INFO.getNumberOfTaxa() + "").
                    replaceFirst("@", PER_TAXA_INFO.getNumberOfBipartitions() + "");
            out.write(header);
        } else if (PHYLIP.equalsIgnoreCase(FORMAT)) {
            String header = PHYLIP_HEADER.replaceFirst("@", PER_TAXA_INFO.getNumberOfTaxa() + "").
                    replaceFirst("@", PER_TAXA_INFO.getNumberOfBipartitions() + "");
            out.write(header);
        }
    }

    private void writeSequenceNameToFile(BufferedWriter out) throws IOException {
        if (NEXUS.equals(FORMAT)) {
            out.write("\t'" + PER_TAXA_INFO.getCurrentSeqName() + "' ");
        } else if (FORMAT.equalsIgnoreCase(PHYLIP)) {
            out.write(PER_TAXA_INFO.getCurrentSeqName() + " ");
        } else {
            out.write(">" + PER_TAXA_INFO.getCurrentSeqName());
            out.newLine();
        }
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: <treesfile> <output> <ouputformat> [-dna] [-randomize seed]\n"
                    + "		<treesfile>: A file containing newick trees, one tree per line\n"
                    + "		<Output>: The name of the output MRP Matrix file\n"
                    + "		<outformat>: use NEXUS for nexus, PHYLIP for phylip, or FASTA for fasta fromatted otuput\n"
                    + " 		-dna: output As and Ts instead of 0 and 1\n"
                    + "		-randomize: randomize 0-1 codings. Seed number is optional.");

            System.exit(1);
        }
        FastMRP mrpCon = new FastMRP(args[0], args[1], args[2]);

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-dna")) {
                mrpCon.setCharacters('A', 'T', '-');
            } else if (args[i].equals("-randomize")) {
                try {
                    int seed = Integer.parseInt(args[i + 1]);
                    mrpCon.random = new Random(seed);
                } catch (NumberFormatException e) {
                    mrpCon.random = new Random();
                } catch (IndexOutOfBoundsException e) {
                    mrpCon.random = new Random();
                }
            }
        }
        try {
            mrpCon.readTreesFile();
            mrpCon.writeMRPToFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private class InfoPerTaxa {

        Integer currentSeqIndex = -1;
        Integer currentColumnInd = -1;
        /*
         * Index of 1st List ~ sequences, 2nd List's values ~ column indecies
         */
        ArrayList<ArrayList<Integer>> columnsPerSequence = new ArrayList<>();
        /*
         * Index of Arraylist ~ sequences, Elements in HashSet ~ Trees
         */
        // TODO: can make this more memory efficient by turning the set to a list
        ArrayList<HashSet<Integer>> treesPerSequence = new ArrayList<>();
        /*
         * Mapping taxa names to IDs 
         */
        HashMap<String, Integer> taxaNameToIndex = new HashMap<>();
        ArrayList<String> taxaNames = new ArrayList<>();
        /*
         * Mapping columns to either 0 or 1 (randomly) to ensure equal 0s or 1s
         * Index ~ 
         */
        private ArrayList<Boolean> columnCoding = new ArrayList<>();

        private Iterator<String> namesIter;
        private Iterator<HashSet<Integer>> treesIter;
        private Iterator<ArrayList<Integer>> columnsIter;

        private Iterator<Integer> currentSeqColumnsIter;
        private HashSet<Integer> currentTrees;
        private ArrayList<Integer> currentColumns;
        private String currentSeqName;

        public Integer mapNamesToIds(String taxon) {
            if (taxaNameToIndex.containsKey(taxon)) {
                return taxaNameToIndex.get(taxon);
            } else {
                currentSeqIndex++;
                taxaNameToIndex.put(taxon, currentSeqIndex);
                taxaNames.add(currentSeqIndex, taxon);
                return currentSeqIndex;
            }
        }

        public int addBipartition(Collection<Integer> bipartition) {
            currentColumnInd++;
            // Randomly assign codings to columns if -randomize is provided
            columnCoding.add(currentColumnInd, random == null ? true : random.nextBoolean());

            for (Integer sequence : bipartition) {
                columnsPerSequence.get(sequence).add(currentColumnInd);
            }
            return currentColumnInd;
        }

        void addTreeToSequence(Integer seq, Integer tree) {
            // If the sequence is encountered for the first time, add it to datastructures.
            if (seq >= treesPerSequence.size()) {
                treesPerSequence.add(seq, new HashSet<>());
                columnsPerSequence.add(seq, new ArrayList<>());
            }
            treesPerSequence.get(seq).add(tree);
        }

        void startSeqIteration() {
            treesIter = treesPerSequence.iterator();
            columnsIter = columnsPerSequence.iterator();
            namesIter = taxaNames.iterator();
        }

        boolean hasMoreTaxa() {
            return treesIter.hasNext();
        }

        void nextSequence() {
            currentTrees = treesIter.next();
            currentColumns = columnsIter.next();
            currentSeqColumnsIter = currentColumns.iterator();
            currentSeqName = namesIter.next();
        }

        boolean currentSeqIsInTree(Integer treeId) {
            return currentTrees.contains(treeId);
        }

        public Integer nextBipartitionIndexForCurrentSequence() {
            if (currentSeqColumnsIter.hasNext()) {
                return currentSeqColumnsIter.next();
            } else {
                return -1;
            }
        }

        public String getCurrentSeqName() {
            return currentSeqName;
        }

        public Integer getNumberOfTaxa() {
            return currentSeqIndex + 1;
        }

        public Integer getNumberOfBipartitions() {
            return currentColumnInd + 1;
        }
    }

    private class TreeEndIndex {
        /*
         * Index to list~ tree, Integer ~ the index of the last column of the tree
         */

        ArrayList<Integer> treeEndIndex = new ArrayList<>();
        private Iterator<Integer> iter;

        public void addEndIndex(int lastColumnIndex) {
            treeEndIndex.add(lastColumnIndex);
        }

        void restartTreeIteration() {
            iter = treeEndIndex.iterator();
        }

        int getNexTreeEndInd() {
            return iter.next();
        }
    }
}
