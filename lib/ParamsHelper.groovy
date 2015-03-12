class ParamsHelper {
    /**
     * Check that the given value is non-null and not the empty string.
     *
     * If this is the case then print a message and exit.
     *
     * @param value Value to be checked for being empty.
     * @param name  Name of the parameter to check.
     */
    static def checkNonEmptyParam(value, name) {
        if (value != null && !"".equals(value))
            return;
        println String.format("Parameter %s empty (%s).", name, value);
        println String.format("Use --%s <value> for setting this parameter.", name);
        System.exit(1);
    }
}
