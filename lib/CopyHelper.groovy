/**
 * Helper class for copying files from channels.
 */
class CopyHelper {
    /** base directory to use for copying */
    public final String baseDir;
    /** whether or not to print copy messages*/
    public final boolean printCopyMsgs;

    /**
     * Initialize object.
     *
     * @param baseDir base directory for copying to.
     * @param printCopyMsgs whether or not to print messages on file copying
     */
    CopyHelper(baseDir, printCopyMsgs) {
        this.baseDir = baseDir;
        this.printCopyMsgs = printCopyMsgs;
    }

    /**
     * Copy files from flattened <code>channel</code> to <code>outPath</code>
     *
     * @param channel source Channel with files, possibles in lists
     * @param outPath output directory, relative to <code>this.baseDir</code>
     * @param unzipZips boolean indicating whether to unzip ".zip" files.
     */
    def copyFiles(channel, outPath, unzipZips = false) {
        def outDir = new File(this.baseDir + "/" + outPath);
        if (!outDir.exists()) {
            if (printCopyMsgs)
                println "mkdir " + outDir;
            outDir.mkdirs();
        }
        channel.flatten().subscribe { f ->
            if (printCopyMsgs)
                println "copy " + f.getName() + " => " + outDir + "/" + f.getName();
            f.copyTo(outDir.toString());
            // maybe unzip file
            if (unzipZips && f.getName().endsWith(".zip")) {
                if (printCopyMsgs)
                    println "unzipping ${f.getName()}";
                Process p = "unzip -d ${outDir} ${outDir}/${f.getName()}".execute();
            }
        }
    }
}
