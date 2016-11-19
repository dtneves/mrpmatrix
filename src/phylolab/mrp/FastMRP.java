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
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
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
    private final char ONE;
    private final char ZERO;
    private final char MISSING;
    /**
     * First list is the STACK position, each Collection is a bipartition (i.e
     * index of taxa in ONE of the partitions)
     */
    private final LinkedList<Collection<Integer>> STACK;
    private final InfoPerTaxa PER_TAXA_INFO;
    private final TreeEndIndex TREE_END_INDEX;
    private Random random;
    private int treeCount;

    public FastMRP(
            final String IN_FILENAME,
            final String OUT_FILENAME,
            final String FORMAT,
            final boolean DNA,
            final boolean RANDOMIZE,
            final Long SEED) {
        this.TREES_FILENAME = IN_FILENAME;
        this.MRP_FILENAME = OUT_FILENAME;
        this.FORMAT = FORMAT;
        if (DNA) {
            this.ONE = 'A';
            this.ZERO = 'T';
            this.MISSING = '-';
        } else {
            this.ONE = '1';
            this.ZERO = '0';
            this.MISSING = '?';
        }
        if (RANDOMIZE) {
            if (SEED == null) {
                this.random = new Random();
            } else {
                this.random = new Random(SEED);
            }
        }
        this.STACK = new LinkedList<>();
        this.PER_TAXA_INFO = new InfoPerTaxa();
        this.TREE_END_INDEX = new TreeEndIndex();
        this.treeCount = 0;
    }

    private void readTreesFile() throws IOException {
        try (BufferedReader in = new BufferedReader(
                new InputStreamReader(new FileInputStream(TREES_FILENAME)))) {
            String tree;
            int lastColumnIndex = -1;

            while ((tree = in.readLine()) != null) {
                final NewickTokenizer TOKENIZER
                        = new NewickTokenizer(tree.trim());

                if (!"(".equals(TOKENIZER.nextToken())) {
                    throw new RuntimeException(
                            "The tree does not start with a (");
                }
                while (TOKENIZER.hasNext()) {
                    String token = TOKENIZER.nextToken();

                    if (token != null) {
                        switch (token) {
                            case "(":
                                STACK.addLast(new ArrayList<>());
                                break;
                            case ")":
                                if (STACK.size() > 0) {
                                    Collection<Integer> top = STACK.getLast();
                                    lastColumnIndex
                                            = PER_TAXA_INFO.addBipartition(top);
                                    STACK.removeLast();
                                    if (STACK.size() > 0) {
                                        STACK.getLast().addAll(top);
                                    }
                                }
                                break;
                            case ";":
                                TREE_END_INDEX.addEndIndex(lastColumnIndex);
                                break;
                            default:
                                Integer seqId
                                        = PER_TAXA_INFO.mapNamesToIds(token);
                                if (STACK.size() > 0) {
                                    STACK.getLast().add(seqId);
                                }
                                PER_TAXA_INFO.addTreeToSequence(
                                        seqId, this.treeCount);
                                break;
                        }
                    }
                }
                this.treeCount++;
            }
        }
    }

    private void writeMRPToFile() throws IOException {
        try (BufferedWriter OUT = new BufferedWriter(
                new FileWriter(MRP_FILENAME))) {
            writeHeaderToFile(OUT);
            // iterate over sequences, ONE by ONE
            PER_TAXA_INFO.startSeqIteration();
            while (PER_TAXA_INFO.hasMoreTaxa()) {
                Integer nextOneColumn;

                PER_TAXA_INFO.nextSequence();
                writeSequenceNameToFile(OUT);
                nextOneColumn = PER_TAXA_INFO
                        .nextBipartitionIndexForCurrentSequence();
                TREE_END_INDEX.restartTreeIteration();
                for (int tree = 0, column = 0; tree < treeCount; tree++) {
                    int treeEndIndex = TREE_END_INDEX.getNexTreeEndInd();

                    if (PER_TAXA_INFO.currentSeqIsInTree(tree)) {
                        while (column <= treeEndIndex) {
                            Boolean coding
                                    = PER_TAXA_INFO.COLUMN_CODING.get(column);

                            if (nextOneColumn == column) {
                                OUT.write(coding ? ONE : ZERO);
                                nextOneColumn = PER_TAXA_INFO
                                        .nextBipartitionIndexForCurrentSequence();
                            } else {
                                OUT.write(coding ? ZERO : ONE);
                            }
                            column++;
                        }
                    } else {
                        char[] missings = new char[treeEndIndex - column + 1];

                        Arrays.fill(missings, MISSING);
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

    private void writeFooterToFile(final BufferedWriter OUT)
            throws IOException {
        if (NEXUS.equalsIgnoreCase(FORMAT)) {
            OUT.write(NEXUS_FOOTER);
        }
    }

    private void writeHeaderToFile(final BufferedWriter OUT)
            throws IOException {
        if (NEXUS.equalsIgnoreCase(FORMAT)) {
            OUT.write(NEXUS_HEADER.replaceFirst("@", PER_TAXA_INFO.getNumberOfTaxa() + "").
                    replaceFirst("@", PER_TAXA_INFO.getNumberOfBipartitions() + ""));
        } else if (PHYLIP.equalsIgnoreCase(FORMAT)) {
            OUT.write(PHYLIP_HEADER.replaceFirst("@", PER_TAXA_INFO.getNumberOfTaxa() + "").
                    replaceFirst("@", PER_TAXA_INFO.getNumberOfBipartitions() + ""));
        }
    }

    private void writeSequenceNameToFile(final BufferedWriter OUT)
            throws IOException {
        if (FORMAT.equals(NEXUS)) {
            OUT.write("\t'" + PER_TAXA_INFO.getCurrentSeqName() + "' ");
        } else if (FORMAT.equalsIgnoreCase(PHYLIP)) {
            OUT.write(PER_TAXA_INFO.getCurrentSeqName() + " ");
        } else {
            OUT.write(">" + PER_TAXA_INFO.getCurrentSeqName());
            OUT.newLine();
        }
    }

    private static void usage() {
        System.err.println("Usage: <trees_file> <output> <output_format> [-dna] [-randomize seed]\n"
                + "		<trees_file>: A file containing Newick trees, one tree per line\n"
                + "		<output>: The name of the output matrix representation (MR) file\n"
                + "		<output_format>: use NEXUS for nexus, PHYLIP for phylip, or FASTA for fasta formatted output\n"
                + " 	-dna: output As and Ts instead of 0 and 1\n"
                + "		-randomize: randomize 0-1 codings, the seed number is optional.");
    }

    public static void main(final String[] ARGS) {
        if (ARGS.length < 3) {
            usage();
        } else {
            final List<String> ARGUMENTS = new ArrayList<>();
            boolean dna = false;
            boolean randomize = false;
            Long seed = null;

            for (int i = 0; i < ARGS.length; i++) {
                if (ARGS[i].equals("-dna")) {
                    dna = true;
                } else if (ARGS[i].equals("-randomize")) {
                    randomize = true;
                    if (i < ARGS.length - 1) {
                        try {
                            seed = Long.parseLong(ARGS[++i]);
                        } catch (final NumberFormatException NFE) {
                            System.err.println(
                                    "The given seed is NOT an integer number.");
                            System.exit(-1);
                        }
                    }
                } else {
                    ARGUMENTS.add(ARGS[i]);
                }
            }
            if (ARGUMENTS.size() < 3) {
                usage();
            } else {
                final FastMRP FAST_MRP = new FastMRP(
                        ARGUMENTS.get(0), ARGUMENTS.get(1), ARGUMENTS.get(2),
                        dna, randomize, seed);
                
                try {
                    FAST_MRP.readTreesFile();
                    FAST_MRP.writeMRPToFile();
                } catch (final IOException IOE) {
                    System.err.println(IOE);
                }
            }
        }
    }

    private final class InfoPerTaxa {
        
        /**
         * Index of 1st List ~ sequences, 2nd List's values ~ column indices
         */
        private final List<List<Integer>> COLUMNS_PER_SEQUENCE;
        /**
         * Index of List ~ sequences, Elements in Set ~ Trees
         */
        // TODO: can make this more memory efficient by turning the set to a list
        private final List<Set<Integer>> TREES_PER_SEQUENCE;
        /**
         * Mapping taxa names to IDs
         */
        private final Map<String, Integer> TAXA_NAME_TO_INDEX;
        private final List<String> TAXA_NAMES;
        /**
         * Mapping columns to either 0 or 1 (randomly) to ensure equal 0s or 1s
         * Index ~
         */
        private final List<Boolean> COLUMN_CODING;

        private Integer currentSeqIndex;
        private Integer currentColumnInd;

        private Iterator<String> namesIter;
        private Iterator<Set<Integer>> treesIter;
        private Iterator<List<Integer>> columnsIter;

        private Iterator<Integer> currentSeqColumnsIter;
        private Set<Integer> currentTrees;
        private List<Integer> currentColumns;
        private String currentSeqName;

        private InfoPerTaxa() {
            COLUMNS_PER_SEQUENCE = new ArrayList<>();
            TREES_PER_SEQUENCE = new ArrayList<>();
            TAXA_NAME_TO_INDEX = new HashMap<>();
            TAXA_NAMES = new ArrayList<>();
            COLUMN_CODING = new ArrayList<>();
            currentSeqIndex = -1;
            currentColumnInd = -1;
        }
        
        public Integer mapNamesToIds(final String TAXON) {
            if (TAXA_NAME_TO_INDEX.containsKey(TAXON)) {
                return TAXA_NAME_TO_INDEX.get(TAXON);
            } else {
                currentSeqIndex++;
                TAXA_NAME_TO_INDEX.put(TAXON, currentSeqIndex);
                TAXA_NAMES.add(currentSeqIndex, TAXON);
                return currentSeqIndex;
            }
        }

        private int addBipartition(final Collection<Integer> BIPARTITION) {
            currentColumnInd++;
            // Randomly assign codings to columns if -randomize is provided
            COLUMN_CODING.add(
                    currentColumnInd,
                    random == null ? true : random.nextBoolean());

            for (Integer sequence : BIPARTITION) {
                COLUMNS_PER_SEQUENCE.get(sequence).add(currentColumnInd);
            }
            return currentColumnInd;
        }

        private void addTreeToSequence(final Integer SEQ, final Integer TREE) {
            // If the sequence is encountered for the first time, 
            // add it to datastructures.
            if (SEQ >= TREES_PER_SEQUENCE.size()) {
                TREES_PER_SEQUENCE.add(SEQ, new HashSet<>());
                COLUMNS_PER_SEQUENCE.add(SEQ, new ArrayList<>());
            }
            TREES_PER_SEQUENCE.get(SEQ).add(TREE);
        }

        private void startSeqIteration() {
            treesIter = TREES_PER_SEQUENCE.iterator();
            columnsIter = COLUMNS_PER_SEQUENCE.iterator();
            namesIter = TAXA_NAMES.iterator();
        }

        private boolean hasMoreTaxa() {
            return treesIter.hasNext();
        }

        private void nextSequence() {
            currentTrees = treesIter.next();
            currentColumns = columnsIter.next();
            currentSeqColumnsIter = currentColumns.iterator();
            currentSeqName = namesIter.next();
        }

        private boolean currentSeqIsInTree(final Integer TREE_ID) {
            return currentTrees.contains(TREE_ID);
        }

        private Integer nextBipartitionIndexForCurrentSequence() {
            return currentSeqColumnsIter.hasNext()
                    ? currentSeqColumnsIter.next() : -1;
        }

        private String getCurrentSeqName() {
            return currentSeqName;
        }

        private Integer getNumberOfTaxa() {
            return currentSeqIndex + 1;
        }

        private Integer getNumberOfBipartitions() {
            return currentColumnInd + 1;
        }
    }

    private final class TreeEndIndex {

        /**
         * Index to list~ tree, Integer ~ the index of the last column of the
         * tree
         */
        private final List<Integer> TREE_END_INDEX;
        private Iterator<Integer> iter;

        private TreeEndIndex() {
            TREE_END_INDEX = new ArrayList<>();
        }

        private void addEndIndex(int lastColumnIndex) {
            TREE_END_INDEX.add(lastColumnIndex);
        }

        private void restartTreeIteration() {
            iter = TREE_END_INDEX.iterator();
        }

        private int getNexTreeEndInd() {
            return iter.next();
        }
    }
}
