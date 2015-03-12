class PathUtil {
	/**
	 * @param left path to left file name
	 * @param right path to right file name
	 * @return longest common prefix of <code>left</code> and <code>right</code>
	 */
	static def lcp(left, right) {
		int end = Math.max(left.size(), right.size());
		int i = 0;
		for (; i < end; ++i)
			if (left[i] != right[i])
				break;
		return left.substring(0, i);
	}
}