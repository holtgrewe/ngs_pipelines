class ChannelUtil {
	/**
	 * Join the files from two Channels and add a label.
	 *
	 * @param label String with a label for the pair, first element of triple.
	 * @param left  left Channel
	 * @param right right Channel
	 * @return new Channel of triplet lists [ label, leftFile, rightFile ]
	 */
	static def createFilePairChannel(label, left, right) {
		return left.merge(right) { l, r -> [ label, l, r ] }
	}
}
