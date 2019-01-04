import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import java.io.*;
import java.util.Dictionary;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;

public class Annotation {

    double maf_threshold, conf_threshold;
    int col_INFO, passQC, passVariants, total_population_sample, total_population_variants, total_vcf_sample, total_vcf_variants, failed_variants;
    String population_file, vcf, tag;
    public Annotation(String argv0, String argv1, String argv2, String argv3) {
        this.maf_threshold = Double.parseDouble(argv0);
        this.conf_threshold = Double.parseDouble(argv1);
        this.population_file = argv2;
        this.vcf = argv3;
        this.passQC = 0;
        this.passVariants = 0;
        this.total_population_sample = 0;
        this.total_population_variants = 0;
        this.total_vcf_sample = 0;
        this.total_vcf_variants = 0;
        this.failed_variants = 0;
        this.tag = "";
    }

    public static class Args {
        @Parameter(names = { "--allele_frequency", "-af" }, description = "Keep variants with allele frequency <= threshold", required = true)
        private String af;

        @Parameter(names = { "--confidence", "-conf" }, description = "Keep variants with confidence > threshold", required = true)
        private String conf;

        @Parameter(names = { "--population", "-p" }, description = "Input population gvcf file", required = true)
        private String p;

        @Parameter(names = { "--vcf", "-v" }, description = "Input vcf file to be filtered", required = true)
        private String v;

        @Parameter(names = { "--export_population", "-e" }, description = "Export population gvcf file with GTF(genotype frequency), AF(allele frequency) and CONF(confidence) annotation", required = false)
        private boolean e;
    }

    private static String passMAFLines(Annotation annotation, String line, String passLine) {
        String out_line = "";
        String[] splits = line.split("\t");
        String[] alt = splits[4].split(",");
        String[] pass_tag = passLine.split("-");
        String passAlt = pass_tag[3];
        Integer passAltID = 0;
        String tag = "";
        for (int i = 0; i < pass_tag.length - 1; i++)
            tag += pass_tag[i] + "-";

        if (!annotation.tag.equals(tag))
            annotation.passVariants += 1;

        for (int i = 0; i < alt.length; i++)
            if (passAlt.equals(alt[i])) passAltID = i;
        out_line = line + "\t" + passAltID;
        annotation.tag = tag;
//        System.out.println(out_line);
        return out_line;
    }

    private static String filteringMAFformat(String line) {
        String out_line = "";
        String[] splits = line.split("\t");
        String[] alt = splits[4].split(",");
        String chr = splits[0];
        String pos = splits[1];
        String ref = splits[3];

        for (int i = 0; i < alt.length; i++) {
            out_line += chr + "-" + pos + "-" + ref + "-" + alt[i] + "\n";
        }
        return out_line.trim();
    }

    // variants maf > threshold
    private static String filteringMAF(double maf, double conf_threshold, String line) {
        String out_line = "";
        String[] splits = line.split("\t");
        String[] alt = splits[4].split(",");
        String chr = splits[0];
        String pos = splits[1];
        String ref = splits[3];
        String conf = splits[splits.length-3];

        if (splits[splits.length-1].equals("-")) return "empty";
        String[] AF = splits[splits.length-1].split(";");
        int tmp = 0;

        for (int i = 0; i < AF.length; i++) {
            String[] text = AF[i].split(":");
            if (Double.parseDouble(text[1]) > maf && Double.parseDouble(conf) > conf_threshold) {
                tmp = Integer.parseInt(text[0]) - 1;
                out_line += chr + "\t" + pos + "\t" + ref + "\t" + alt[tmp] + "\n";
            }
        }
        if (out_line.equals("")) return "empty";
        return out_line.trim();
    }

    // count AF from GTF
    private static String AlleleFrequencyCounter(String line) {
        String out_line = "";
        String[] splits = line.split("\t");
        String[] GTfreq = splits[splits.length-1].split(";");
        Dictionary AF = new Hashtable();
        String[] allele;
//        System.out.println(splits[splits.length-1]);

        for (int i = 0; i < GTfreq.length; i++) {
            String GT = GTfreq[i].split(":")[0];
            if (GT.contains("|"))
                allele = GT.split("\\|");
            else
                allele = GT.split("/");

            try {
                for (int j = 0; j < allele.length; j++) {
                    if (!allele[j].equals("0") && !((Hashtable) AF).containsKey(Integer.parseInt(allele[j])))
                        AF.put(Integer.parseInt(allele[j]), 0.5 * Double.parseDouble(GTfreq[i].split(":")[1]));

                    else if (!allele[j].equals("0") && ((Hashtable) AF).containsKey(Integer.parseInt(allele[j])))
                        AF.put(Integer.parseInt(allele[j]), (double) AF.get(Integer.parseInt(allele[j])) + (0.5 * Double.parseDouble(GTfreq[i].split(":")[1])));
                }
            } catch (NumberFormatException nfe) {
                out_line = "-;";
                return line + "\t" +out_line.replaceFirst(".$","");
            }
        }

        Iterator it = ((Hashtable) AF).entrySet().iterator();
        while(it.hasNext()){
            Map.Entry me = (Map.Entry)it.next();
            out_line += me.getKey() + ":" + me.getValue() + ";";
        }

        if (out_line.equals("")) out_line = "-;";
        return line + "\t" +out_line.replaceFirst(".$","");
    }

    private static String GenotypeFrequencyCounter(int num_info, int num_sample, String line) {
        String out_line = "";
        String[] splits = line.split("\t");
        Dictionary GTF = new Hashtable();
        double num_hasINFO = 0.0;
        int tmp = num_sample - num_info;

        if (num_sample != splits.length - num_info) return "failed";

        for (int i = splits.length - num_sample; i < splits.length; i++) {
            String GT = splits[i].split(":")[0];

//            System.out.print(GT+"\t");
            // drop GT contains "." (without info.)
            if (!((Hashtable) GTF).containsKey(GT) && !GT.contains(".")) {
                GTF.put(GT, 1.0);
                num_hasINFO += 1.0;
            }
            else if (((Hashtable) GTF).containsKey(GT) && !GT.contains(".")){
                GTF.put(GT, (double) GTF.get(GT) + 1.0);
                num_hasINFO += 1.0;
            }
        }

        Iterator it = ((Hashtable) GTF).entrySet().iterator();
        while(it.hasNext()){
            Map.Entry me = (Map.Entry)it.next();
            out_line += me.getKey() + ":" + (double) me.getValue() / num_hasINFO + ";";
        }
//        System.out.println(out_line);
//        System.out.println(num_hasINFO / (double) num_sample + "\t" + out_line.replaceFirst(".$",""));
        return line + "\t" + num_hasINFO / (double) num_sample + "\t" +out_line.replaceFirst(".$","");
    }


    private static int[] sampleCounter(String line) {
        int[] ct_sample = new int[2];
        ct_sample[0] = 0;
        ct_sample[1] = 0;
        String[] samples = line.split("\t");
        for (int i = 0; i < samples.length; i++) {
            ct_sample[0]++;
            if (samples[i].endsWith("FORMAT")) {
                ct_sample[1] = ct_sample[0];
                ct_sample[0] = samples.length - ct_sample[1];
                return ct_sample;
            }
        }
        return ct_sample;
    }


    public static void main(String[] args) throws IOException {

        Args argv = new Args();
        new JCommander(argv).parse(args);

        long startTime = System.currentTimeMillis();

        Annotation annotation = new Annotation(argv.af, argv.conf, argv.p, argv.v);
        Dictionary mapper = new Hashtable();

        // population gvcf
        try (BufferedReader br = new BufferedReader(new FileReader(annotation.population_file))) {
            String line, line_GTF, line_AF, line_filtered;
            String[] line_filtered_array;

            FileWriter fileWriter = null;
            PrintWriter printWriter = null;
            if (argv.e) {
                fileWriter = new FileWriter(annotation.population_file.split("\\.")[0] + "_GTF_AF.g.vcf");
                printWriter = new PrintWriter(fileWriter);
            }
            while ((line = br.readLine())!=null) {
                if (annotation.total_population_variants % 10000 == 0 && annotation.total_population_variants != 0)
                    System.out.println("loading population database: " + annotation.total_population_variants + " variants");

                if (line.startsWith("#")) {
                    int[] tmp = sampleCounter(line);
                    annotation.total_population_sample = tmp[0];
                    annotation.col_INFO = tmp[1];
                    if (argv.e)
                        printWriter.println(line+"\tCONF"+"\tGTF"+"\tAF");
                }
                else {
//                    System.out.println(line);
                    line_GTF = GenotypeFrequencyCounter(annotation.col_INFO, annotation.total_population_sample, line);
                    if (!line_GTF.equals("failed")) {
                        line_AF = AlleleFrequencyCounter(line_GTF);
                        if (argv.e)
                            printWriter.println(line_AF);
                        line_filtered = filteringMAF(annotation.maf_threshold, annotation.conf_threshold, line_AF);
                        if (line_filtered != "empty") {
                            line_filtered_array = line_filtered.split("\n");
                            for (int i = 0; i < line_filtered_array.length; i++) {
                                mapper.put(line_filtered_array[i].replace("\t", "-"), 1);
                            }
                        }
                    }
                    else
                        annotation.failed_variants++;
                    annotation.total_population_variants++;
                }
            }
            int tmp = annotation.total_population_variants - annotation.failed_variants;
            System.out.println("INFO\t[ population database ] has: " + annotation.total_population_sample + " samples");
            System.out.println("INFO\ttotal population variants: " + annotation.total_population_variants + ", with " + tmp + " loaded variants, and " + annotation.failed_variants + " failed variants");
            if (argv.e)
                printWriter.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // vcf
        try (BufferedReader br = new BufferedReader(new FileReader(annotation.vcf))) {
            String line, line_filtered;
            String[] line_filtered_array;

            FileWriter fileWriter = new FileWriter(annotation.vcf.split("\\.")[0] + "_filtered_maf_below_" + annotation.maf_threshold + "_conf_above_" + annotation.conf_threshold + ".vcf");
            PrintWriter printWriter = new PrintWriter(fileWriter);
            while ((line = br.readLine())!=null) {
                if (annotation.total_vcf_variants % 5000 == 0 && annotation.total_vcf_variants != 0)
                    System.out.println("processing vcf: " + annotation.total_vcf_variants + " variants");

                if (line.startsWith("#")) {
                    annotation.total_vcf_sample = sampleCounter(line)[0];
                    printWriter.println(line+"\tfiltered_maf_above_" + annotation.maf_threshold);
                }
                else {
                    line_filtered = filteringMAFformat(line);
                    line_filtered_array = line_filtered.split("\n");
                    for (int i = 0; i < line_filtered_array.length; i++) {
                        if (!((Hashtable) mapper).containsKey(line_filtered_array[i]))
                            printWriter.println(passMAFLines(annotation, line, line_filtered_array[i]));
                    }
                    annotation.total_vcf_variants++;
                }
            }
            printWriter.close();

            FileWriter logfileWriter = new FileWriter(annotation.vcf.split("\\.")[0] + "_filtered_maf_below_" + annotation.maf_threshold + "_conf_above_" + annotation.conf_threshold + ".log");
            PrintWriter printlogWriter = new PrintWriter(logfileWriter);

            int tmp = annotation.total_population_variants - annotation.failed_variants;
            System.out.println("INFO\tallele frequency threshold: " + argv.af + ", confidence threshold: " + argv.conf);
            System.out.println("INFO\t[ population database ] has: " + annotation.total_population_sample + " samples");
            System.out.println("INFO\ttotal population variants: " + annotation.total_population_variants + ", with " + tmp + " loaded variants, and " + annotation.failed_variants + " failed variants");
            System.out.println("INFO\t[ vcf ] has: " + annotation.total_vcf_sample + " samples");
            System.out.println("INFO\tpass/total variants: " + annotation.passVariants + "/" + annotation.total_vcf_variants);
            System.out.println("INFO\tUsing Time: " + (System.currentTimeMillis() - startTime)/1000 + " seconds");

            printlogWriter.println("INFO\tallele frequency threshold: " + argv.af + ", confidence threshold: " + argv.conf);
            printlogWriter.println("INFO\t[ population database ] has: " + annotation.total_population_sample + " samples");
            printlogWriter.println("INFO\ttotal population variants: " + annotation.total_population_variants + ", with " + tmp + " loaded variants, and " + annotation.failed_variants + " failed variants");
            printlogWriter.println("INFO\t[ vcf ] has: " + annotation.total_vcf_sample + " samples");
            printlogWriter.println("INFO\tpass/total variants: " + annotation.passVariants + "/" + annotation.total_vcf_variants);
            printlogWriter.println("INFO\tUsing Time: " + (System.currentTimeMillis() - startTime)/1000 + " seconds");
            printlogWriter.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
