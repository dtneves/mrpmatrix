package phylolab;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class NewickTokenizer {

    private final Matcher PATTERN;
    private final boolean STRIP;

    public NewickTokenizer(final String INPUT) {
        this.STRIP = true;
        if (this.STRIP) {
            PATTERN = Pattern.compile(
                    "([(])|([)][^,:;)]*)|([;])|(:)|([^,);(:]*)").matcher(INPUT);
        } else {
            PATTERN = Pattern.compile(
                    "([(])|([)][^,;)]*)|([;])|([^,);(]*)").matcher(INPUT);
        }
        PATTERN.find();
    }

    public NewickTokenizer(final String INPUT, final boolean STRIP) {
        this.STRIP = STRIP;
        if (this.STRIP) {
            PATTERN = Pattern.compile(
                    "([(])|([)][^,:;)]*)|([;])|(:)|([^,);(:]*)").matcher(INPUT);
        } else {
            PATTERN = Pattern.compile(
                    "([(])|([)][^,;)]*)|([;])|([^,);(]*)").matcher(INPUT);
        }
        PATTERN.find();
    }

    public final boolean hasNext() {
        return !PATTERN.hitEnd();
    }

    public final String nextToken() {
        final String RES = PATTERN.group();

        PATTERN.find();
        // This is to STRIP off any support value / internal label nodes 
        // that can follow a left bracket.
        if (STRIP && RES.startsWith(")")) {
            return ")";
        }
        if (RES != null) {
            switch (RES) {
                case "":
                    final String TOKEN = nextToken();
                    
                    return TOKEN.isEmpty() ? null : TOKEN;
                case ":":
                    // This is to STRIP off the branch lenght values.
                    nextToken();
                    return nextToken();
            }
        }
        return RES != null ? RES.trim() : null;
    }

    public static void main(final String[] ARGS) {
        //String tree = "(t90:0.01,((t18:0.04,(((t17:0.0,t49:0.09):0.00,(((((t81:0.037,t45:0.02):0.02,t66:0.01):0.00724936580827239413,t41:0.05336582750396093311):0.03045292104932446550,(t14:0.11614645744250107207,((t16:0.05973131223865356387,t91:0.06115882420120943158):0.03367781109364907655,((t30:0.06168112199284891961,(t82:0.00707911640722609301,t19:0.01702795908063739483):0.04105036384987210962):0.04367847525636570777,(t88:0.08520512585491028801,((t77:0.02231272403400484661,t100:0.03396176535082153641):0.03970879193218978392,t46:0.07474947196269035588):0.03491872410661160664):0.00983925300854607623):0.00000122920446704147):0.00271415429521279158):0.00799815476646543837):0.01024163042108295132,(t80:0.07155782866412271903,t68:0.09571882482752180898):0.01525096682579810646):0.01622216107791855585):0.00574275591113855306,t31:0.08260865518257545781):0.00767107143981279709):0.00100486256845222326,t15:0.08070670351238432016):0.03897657131187436119,t43:0.01734246018246490828):0.0;";
        String tree = "(2sd:32[2],((M:29[23],(A:10[2],(C:1,D:3)0.30[12]:232,B:12),N)22:232[12],L));";
        NewickTokenizer tokenizer = new NewickTokenizer(tree, false);

        while (tokenizer.hasNext()) {
            System.out.print(tokenizer.nextToken() + "  ");
        }
    }

}
