package cmsc420_s23;

import java.util.*;

public class SMkdTree<LPoint extends LabeledPoint2D> {




	private abstract class Node { // generic node type

		InternalNode parent; // our parent
		LinkedList<LPoint> contenders;

		Node(){
			this.parent = null;
			this.contenders = new LinkedList<>();
		}

		// -----------------------------------------------------------------
		// Standard dictionary helpers
		// -----------------------------------------------------------------

		abstract LPoint find(Point2D q); // find point in subtree

		abstract Node insert(LPoint pt, Rectangle2D cell) throws Exception; // insert a single point

		abstract Node delete(Point2D pt) throws Exception; // delete point

		abstract ArrayList<String> list(); // list entries in reverse preorder

		abstract LPoint nearNeigh(Point2D center, Rectangle2D cell, LPoint best, ArrayList<LPoint> visited); // nearest-neighbor helper

		// -----------------------------------------------------------------
		// Other utilities
		// -----------------------------------------------------------------

		abstract String debugPrint(String prefix); // print for debugging

		abstract public String toString(); // string representation

		abstract ArrayList<LPoint> bruteForceNN(Point2D center); // brute-force nearest neighbor

		abstract Node rebuildAfterInsert(LPoint pt, Rectangle2D cell) throws Exception; // check for rebuild

		abstract void buildPointList(List<LPoint> pts); // build a list of points in this subtree

		abstract Node rebuild(Rectangle2D cell) throws Exception; // rebuild subtree rooted at this node

		abstract int getSize(); // get number of points in subtree

		abstract ExternalNode leftMost(); // get leftmost descendant

		// newly added
		abstract ArrayList<LPoint> addCenter(LPoint center, Rectangle2D cell) throws Exception;
		abstract ArrayList<String> listWithCenters();


		public abstract ArrayList<AssignedPair> listAssignments();

		public abstract LPoint closestCenter(LPoint point);

		//public abstract ArrayList<LPoint> getSites(LabeledPoint2D center, Rectangle2D cell, ArrayList<LPoint> sites);
	}

	// -----------------------------------------------------------------
	// Internal node
	// -----------------------------------------------------------------

	private class InternalNode extends Node {

		int cutDim; // the cutting dimension (0 = x, 1 = y)
		double cutVal; // the cutting value
		Node left, right; // children
		int size; // number of points
		int insertCt; // insert counter

		/**
		 * Constructor
		 */
		InternalNode(int cutDim, double cutVal) {
			this.parent = null;
			this.cutDim = cutDim;
			this.cutVal = cutVal;
			this.left = null;
			this.right = null;
			this.insertCt = 0;
			this.size = 0;
		}

		// -----------------------------------------------------------------
		// Internal-node Utilities
		// -----------------------------------------------------------------

		boolean onLeft(Point2D pt) { // in the left subtree? (for Point2D)
			return pt.get(cutDim) < cutVal;
		}

		boolean onLeft(LPoint pt) { // in the left subtree? (for LPoint)
			return pt.get(cutDim) < cutVal;
		}

		Rectangle2D leftPart(Rectangle2D cell) { // left part of cell
			return cell.leftPart(cutDim, cutVal);
		}

		Rectangle2D rightPart(Rectangle2D cell) { // right part of cell
			return cell.rightPart(cutDim, cutVal);
		}

		// -----------------------------------------------------------------
		// Internal-node helpers
		// -----------------------------------------------------------------

		/**
		 * Find point in this subtree
		 */
		LPoint find(Point2D q) {
			if (onLeft(q)) return left.find(q);
			else return right.find(q);
		}

		/**
		 * Insert a single point.
		 */
		Node insert(LPoint pt, Rectangle2D cell) throws Exception{
			if (onLeft(pt)) {
				left = left.insert(pt, leftPart(cell));
				left.parent = this;
			} else {
				right = right.insert(pt, rightPart(cell));
				right.parent = this;
			}
			size += 1; // increment our size
			insertCt += 1; // increment our insert counter
			return this;
		}

		/**
		 * Check/rebuild if needed after insertion.
		 */
		Node rebuildAfterInsert(LPoint pt, Rectangle2D cell) throws Exception {
			if (insertCt > (size + rebuildOffset)/2) { // time to rebuild
				return rebuild(cell);
			} else if (onLeft(pt)) {
				left = left.rebuildAfterInsert(pt, leftPart(cell));
				left.parent = this;
			} else {
				right = right.rebuildAfterInsert(pt, rightPart(cell));
				right.parent = this;
			}
			return this;
		}

		/**
		 * Delete a point from this subtree. We recurse on the side
		 * containing the point. A null value signifies that our
		 * child is an external node that has lost its last point,
		 * and we respond by returning the other child (since we are
		 * no longer needed).
		 */
		Node delete(Point2D pt) throws Exception {
			if (onLeft(pt)) {
				left = left.delete(pt);
				left.parent = this;
			}
			else {
				right = right.delete(pt);
				right.parent = this;
			}
			size -= 1; // one less point
			return this;
		}

		/**
		 * Get a list of the nodes in reverse preorder.
		 *
		 * Adds current key in parentheses, followed by left and right subtrees.
		 */
		ArrayList<String> list() {
			ArrayList<String> list = new ArrayList<String>();
			list.add(toString()); // add this node
			list.addAll(right.list()); // add right
			list.addAll(left.list()); // add left
			return list;
		}

		/**
		 * Nearest neighbor helper. The closest point seen so far is
		 * best. We search the closer subtree first and the other
		 * subtree if it is still viable.
		 */
		LPoint nearNeigh(Point2D center, Rectangle2D cell, LPoint best, ArrayList<LPoint> visited) {
			// child cells
			Rectangle2D leftCell = cell.leftPart(cutDim, cutVal);
			Rectangle2D rightCell = cell.rightPart(cutDim, cutVal);
			if (onLeft(center)) { // left side is closer
				best = left.nearNeigh(center, leftPart(cell), best, visited);
				double bestDist = (best == null ? Double.POSITIVE_INFINITY : center.distanceSq(best.getPoint2D()));
				if (rightCell.distanceSq(center) <= bestDist) {
					best = right.nearNeigh(center, rightPart(cell), best, visited);
				}
			} else { // right side is closer
				best = right.nearNeigh(center, rightPart(cell), best, visited);
				double bestDist = (best == null ? Double.POSITIVE_INFINITY : center.distanceSq(best.getPoint2D()));
				if (leftCell.distanceSq(center) <= bestDist) {
					best = left.nearNeigh(center, leftPart(cell), best, visited);
				}
			}
			return best;
		}

		/**
		 * Rebuild the subtree rooted at this node. We first compute a sorted list of
		 * external nodes and then rebuild the subtree recursively.
		 */
		Node rebuild(Rectangle2D cell) throws Exception {
			List<LPoint> pts = new ArrayList<LPoint>(); // list for points
			buildPointList(pts); // collect the points in this subtree
			return bulkCreate(pts, cell); // build a new tree
		}

		/**
		 * Builds a list of points in this subtree.
		 */
		void buildPointList(List<LPoint> pts) {
			left.buildPointList(pts); // add left
			right.buildPointList(pts); // add right
		}

		/**
		 * Get leftmost descendant. This is used in the iterator.
		 */
		ExternalNode leftMost() {
			return left.leftMost();
		}


		@Override
		ArrayList<LPoint> addCenter(LPoint center, Rectangle2D cell) throws Exception {
			// check if is in contender list
			ArrayList<LPoint> updated = new ArrayList<>();
			double Rmin = Double. MAX_VALUE;
			for(LPoint curr: this.contenders){
				double currDis = cell.maxDistanceSq(curr.getPoint2D());
				if(currDis< Rmin){
					Rmin = cell.maxDistanceSq(curr.getPoint2D());
				}
			}
//			for(LPoint curr: this.contenders){
//				double Ri =cell.distanceSq(curr.getPoint2D());
//				if( Ri > Rmin){
//					this.contenders.remove(curr);
//				}
//			}

			double dis = cell.distanceSq(center.getPoint2D());
			// if dis less than Rmin add to contenders otherwise not a contender
			if(dis <= Rmin){
				this.contenders.add(center);

				updated.addAll(left.addCenter(center, leftPart(cell)));
				updated.addAll(right.addCenter(center, rightPart(cell)));

			}
			Iterator<LPoint> iterator = this.contenders.iterator();
			for(LPoint curr: this.contenders){
				if(cell.maxDistanceSq(curr.getPoint2D()) < Rmin){
					Rmin = cell.maxDistanceSq(curr.getPoint2D());
				}
			}
			while (iterator.hasNext()) {
				LPoint curr = iterator.next();
				double ri = cell.distanceSq(curr.getPoint2D());
				if (ri > Rmin) {
					// remove current element from contenders list
					iterator.remove();
				}
			}
			return updated;
		}

		@Override
		ArrayList<String> listWithCenters() {
			ArrayList<String> list = new ArrayList<String>();
			Collections.sort(contenders, new LabelComparator());
			list.add(toString() + " => "+getContenders(this.contenders)); // add this node
			list.addAll(right.listWithCenters()); // add right
			list.addAll(left.listWithCenters()); // add left
			return list;
		}

		@Override
		public ArrayList<AssignedPair> listAssignments() {
			ArrayList<AssignedPair> ans = new ArrayList<>();
			ans.addAll(left.listAssignments());
			ans.addAll(right.listAssignments());
			return ans;
		}

		@Override
		public LPoint closestCenter(LPoint q) {
			if (onLeft(q)) return left.closestCenter(q);
			else return right.closestCenter(q);
		}
//
//		@Override
//		public  ArrayList<LPoint> getSites(LabeledPoint2D center, Rectangle2D cell, ArrayList<LPoint> sites) {
//			return null;
//		}



		String getContenders(LinkedList<LPoint> contender){
			StringBuilder sb = new StringBuilder("{");
			int size = contender.size();
			if (size > 10) {
				List<String> firstTenLabels = new ArrayList<>();
				for (int i = 0; i < 10; i++) {
					firstTenLabels.add(contender.get(i).getLabel());
					sb.append(contenders.get(i).getLabel());
					if(i < 9){
						sb.append(" ");
					}
				}
				sb.append("...");
			} else {
				for (int i = 0; i < size; i++) {
					sb.append(contender.get(i).getLabel());
					if (i < size - 1) {
						sb.append(" ");
					}
				}
			}
			sb.append("}");
			return sb.toString();

		}

		// -----------------------------------------------------------------
		// Internal-node utilities
		// -----------------------------------------------------------------

		/**
		 * Get subtree size.
		 */
		int getSize() {
			return size;
		}

		/**
		 * String representation.
		 */
		public String toString() {
			String cutIndic = (cutDim == 0 ? "x=" : "y=");
			return "(" + cutIndic + cutVal + ") " + size + ":" + insertCt;
		}

		/**
		 * Brute-force nearest neighbor (for validation)
		 */
		ArrayList<LPoint> bruteForceNN(Point2D center) {
			ArrayList<LPoint> leftClose = left.bruteForceNN(center); // get left and right closest
			ArrayList<LPoint> rightClose = right.bruteForceNN(center);
			if (leftClose.isEmpty()) { // if either is empty, return the other
				return rightClose;
			} else if (rightClose.isEmpty()) {
				return leftClose;
			} else { // both are nonempty
				double leftDist = center.distanceSq(leftClose.get(0).getPoint2D());
				double rightDist = center.distanceSq(rightClose.get(0).getPoint2D());
				if (leftDist < rightDist) { // left is closer
					return leftClose;
				} else if (rightDist < leftDist) { // right is closer
					return rightClose;
				} else { // same distances
					leftClose.addAll(rightClose); // merge both lists
					return leftClose;
				}
			}
		}

		/**
		 * Debug print subtree.
		 *
		 * This prints the subtree for debugging.
		 */
		String debugPrint(String prefix) {
			return  right.debugPrint(prefix + "| ") + System.lineSeparator() + // right child
					prefix + toString() + System.lineSeparator() + // this node
					left.debugPrint(prefix + "| "); // left child
		}
	}

	// -----------------------------------------------------------------
	// External node
	// -----------------------------------------------------------------

	private class ExternalNode extends Node {

		LPoint point; // the point (null if no point)
		LPoint nearCenter;


		/**
		 * Default constructor.
		 */
		ExternalNode(LPoint point) {
			this.point = point;
			this.parent = null;
		}

		// -----------------------------------------------------------------
		// External-node helpers
		// -----------------------------------------------------------------

		/**
		 * Find point in external node.
		 */
		LPoint find(Point2D q) {
			if (point != null && point.getPoint2D().equals(q)) {
				return point;
			} else {
				return null;
			}
		}

		/**
		 * Insert a single point.
		 */
		Node insert(LPoint pt, Rectangle2D cell) throws Exception {
			if (point == null) {
				point = pt;
				return this;
			}
			else {
				if (point.getPoint2D().equals(pt.getPoint2D())) {
					throw new Exception("Insertion of duplicate point");
				}
				ArrayList<LPoint> pts = new ArrayList<LPoint>();
				pts.add(point);
				pts.add(pt);
				return bulkCreate(pts, cell); // create new tree with both points
			}
		}

		/**
		 * Rebuild if needed after insert. (No effect on externals.)
		 */
		Node rebuildAfterInsert(LPoint pt, Rectangle2D cell) {
			return this;
		}

		/**
		 * Delete a point from this node.
		 */
		Node delete(Point2D pt) throws Exception {
			if (point == null || !point.getPoint2D().equals(pt)) {
				throw new Exception("Deletion of nonexistent point");
			}
			point = null;
			return this;
		}

		/**
		 * List the contents of this node
		 */
		ArrayList<String> list() {
			ArrayList<String> list = new ArrayList<String>();
			list.add(toString()); // add this node
			return list;
		}

		/**
		 * Nearest neighbor. Returns the closest point to center.
		 */
		LPoint nearNeigh(Point2D center, Rectangle2D cell, LPoint best, ArrayList<LPoint> visited) {
			if (visited != null && this.point != null)
				visited.add(this.point);
			return closerOf(center, best, point);
		}

		/**
		 * Rebuild this subtree. (Does nothing for external nodes.)
		 */
		Node rebuild(Rectangle2D cell) {
			return this;
		}

		/**
		 * Builds a list of points in this subtree.
		 */
		void buildPointList(List<LPoint> pts) {
			if (point != null) pts.add(point);
		}

		/**
		 * Get leftmost descendant.
		 */
		ExternalNode leftMost() {
			return this;
		}


		@Override
		ArrayList<LPoint> addCenter(LPoint center, Rectangle2D cell) throws Exception {
			ArrayList<LPoint> updated = new ArrayList<>();
			double Rmin = Double. MAX_VALUE;
			for(LPoint curr: this.contenders){
				if(cell.maxDistanceSq(curr.getPoint2D()) < Rmin){
					Rmin = cell.maxDistanceSq(curr.getPoint2D());
				}
			}
			double dis = cell.distanceSq(center.getPoint2D());
			 if (dis <= Rmin){
			 	this.contenders.add(center);
//			 	if(this.point == null){
//			 		if(nearCenter == null){
//			 			nearCenter = center;
//			 			updated.add(point);
//					}else{
//						LPoint closer = closerOf(this.point.getPoint2D(), center, nearCenter);
//						if(closer.getY() == center.getY() && closer.getX() == center.getX()){
//							nearCenter = center;


//						}
//					}
//				}
			 }

			LPoint nearCenter = closestCenter(center);
			 if(nearCenter != null && nearCenter.getX() == center.getX() && nearCenter.getY() == center.getY()){
				 updated.add(point);
			 }



			Iterator<LPoint> iterator = this.contenders.iterator();
			for(LPoint curr: this.contenders){
				if(cell.maxDistanceSq(curr.getPoint2D()) < Rmin){
					Rmin = cell.maxDistanceSq(curr.getPoint2D());
				}
			}
			while (iterator.hasNext()) {
				LPoint curr = iterator.next();
				double ri = cell.distanceSq(curr.getPoint2D());
				if (ri > Rmin) {
					// remove current element from contenders list
					iterator.remove();
				}
			}

			return updated;
		}

		@Override
		ArrayList<String> listWithCenters() {
			Collections.sort(contenders, new LabelComparator());
			ArrayList<String> ans = new ArrayList<String>();
			ans.add(toString() +" => " + getContenders(this.contenders));
			return ans;
		}



		@Override
		public ArrayList<AssignedPair> listAssignments() {
			ArrayList<AssignedPair> ans = new ArrayList<>();
			if(point == null){
				return ans;
			}
			AssignedPair pair = new AssignedPair();
			double closest = Double.MAX_VALUE;
			LPoint closestCenter = null;
			for(LPoint curr: this.contenders){
				Double currDis = curr.getPoint2D().distanceSq(point.getPoint2D());
				if (closest >= currDis) {
					closestCenter = curr;
					closest = currDis;
				}
			}
			pair.center = closestCenter;
			pair.site = point;
			pair.distanceSq = closestCenter.getPoint2D().distanceSq(point.getPoint2D());
			ans.add(pair);
			return ans;
		}

		@Override
		public LPoint closestCenter(LPoint point) {
			if(this.point == null){
				return null;
			}
			double closest = Double.MAX_VALUE;
			LPoint closestCenter = null;
			for(LPoint curr: this.contenders){
				Double currDis = curr.getPoint2D().distanceSq(this.point.getPoint2D());
				if (closest >= currDis) {
					closestCenter = curr;
					closest = currDis;
				}
			}
			return closestCenter;
		}

//		@Override
//		public  ArrayList<LPoint> getSites(LabeledPoint2D center, Rectangle2D cell, ArrayList<LPoint> sites) {
//			ArrayList<LPoint> updated = new ArrayList<>();
//			double Rmin = Double. MAX_VALUE;
//			for(LPoint curr: this.contenders){
//				if(cell.maxDistanceSq(curr.getPoint2D()) < Rmin){
//					Rmin = cell.maxDistanceSq(curr.getPoint2D());
//				}
//			}
//			LPoint nearCenter = closestCenter(point);
//			if(nearCenter.getX() == center.getX() && nearCenter.getY() == center.getY()){
//
//			}
//			return sites;
//		}



		String getContenders(LinkedList<LPoint> contender){
			StringBuilder sb = new StringBuilder("{");
			int size = contender.size();
			if (size > 10) {
				List<String> firstTenLabels = new ArrayList<>();
				for (int i = 0; i < 10; i++) {
					firstTenLabels.add(contender.get(i).getLabel());
					sb.append(contenders.get(i).getLabel());
					if(i < 9){
						sb.append(" ");
					}
				}
				sb.append("...");
			} else {
				for (int i = 0; i < size; i++) {
					sb.append(contender.get(i).getLabel());
					if (i < size - 1) {
						sb.append(" ");
					}
				}
			}
			sb.append("}");
			return sb.toString();

		}
		// -----------------------------------------------------------------
		// External-node Utilities
		// -----------------------------------------------------------------

		/**
		 * Get subtree size.
		 */
		int getSize() {
			return point == null ? 0 : 1;
		}

		/**
		 * String representation.
		 */
		public String toString() {
			String contents = new String();
			contents += "[" + (point == null ? "null" : point.toString()) + "]";
			return contents;
		}

		/**
		 * Debug print.
		 */
		String debugPrint(String prefix) {
			return prefix + toString();
		}

		/**
		 * Nearest neighbor by brute force (for validation). Returns
		 * the point in this node if it is non-null.
		 */
		ArrayList<LPoint> bruteForceNN(Point2D center) {
			ArrayList<LPoint> result = new ArrayList<LPoint>();
			if (point != null) result.add(point);
			return result;
		}
	}

	// -----------------------------------------------------------------
	// Iterator
	// -----------------------------------------------------------------

	private class LPointIterator implements Iterator<LPoint> {
		private ExternalNode next; // external node with the next point

		/**
		 * Constructor. If tree is empty, the iterator is null. Otherwise
		 * Initialize to first point in leftmost leaf.
		 */
		public LPointIterator() {
			if (nPoints == 0) {
				next = null;
			} else {
				next = root.leftMost();
			}
		}

		/**
		 * More elements remaining?
		 */
		public boolean hasNext() {
			advanceToNonNull();
			return next != null;
		}

		/**
		 * Return current and advance to next.
		 */
		public LPoint next() throws NoSuchElementException {
			advanceToNonNull();
			if (next == null) {
				throw new NoSuchElementException();
			}
			LPoint result = next.point; // final result
			advanceOnce();
			return result;
		}

		/**
		 * Advance to next non-null point.
		 */
		private void advanceToNonNull() {
			while (next != null && next.point == null) {
				advanceOnce();
			}
		}

		/**
		 * Advance to next (which might be null).
		 */
		private void advanceOnce() {
			Node v = next;
			InternalNode u = v.parent; // move along right chain
			while (u != null && u.right == v) {
				v = u;
				u = u.parent;
			}
			if (u == null) { // hit the root?
				next = null;
			} else {
				next = u.right.leftMost(); // leftmost leaf of right child
			}
		}
	}

	// -----------------------------------------------------------------
	// Tree private data
	// -----------------------------------------------------------------

	private final int rebuildOffset; // offset for rebuild counters
	private Rectangle2D rootCell; // the bounding box
	private Node root; // root of the tree
	private int nPoints; // number of points in the tree
	private int deleteCt; // deletion counter

	private LPoint startCenter;
//	private ArrayList<LPoint> centers;


	// -----------------------------------------------------------------
	// Tree public members
	// -----------------------------------------------------------------

	/**
	 * Creates an empty tree.
	 */
	SMkdTree(int rebuildOffset, Rectangle2D rootCell, LPoint startCenter) throws Exception {
		this.rebuildOffset = rebuildOffset;
		this.rootCell = new Rectangle2D(rootCell);
		this.startCenter = startCenter;
		clear();
	}

	/**
	 * Clear the tree.
	 */
	public void clear() {
		root = new ExternalNode(null);
		nPoints = 0;
		deleteCt = 0;
	}

	/**
	 * Number of points.
	 */
	public int size() {
		return nPoints;
	}

	/**
	 * Delete counter value.
	 */
	public int deleteCount() {
		return deleteCt;
	}

	/**
	 * Find an point in the tree.
	 */
	public LPoint find(Point2D q) {
		return root.find(q);
	}

	/**
	 * Insert a point.
	 */
	public void insert(LPoint pt) throws Exception {
		if (!rootCell.contains(pt.getPoint2D())) { // outside bounds?
			throw new Exception("Attempt to insert a point outside bounding box");
		}
		root = root.insert(pt, rootCell); // insert
		root = root.rebuildAfterInsert(pt, rootCell); // rebuild if needed
		root.parent = null;
		nPoints += 1; // update point count
	}

	/**
	 * Delete a point. Note that the point being deleted does not need to match
	 * fully. It suffices that it has enough information to satisfy the comparator.
	 *
	 * We first check whether the point is in the tree, using find. This is
	 * necessary, because (due to fact that points lying on the splitting line
	 * be in either subtree) it is difficult for the deletion helper to check
	 * for nonexistent points.
	 */
	public void delete(Point2D pt) throws Exception {
		root = root.delete(pt);
		root.parent = null;
		nPoints -= 1; // one fewer point
		deleteCt += 1; // increment deletion counter
		if (deleteCt > nPoints) {
			root = root.rebuild(rootCell); // time to rebuild
			deleteCt = 0; // reset the delete counter
		}
	}

	/**
	 * Get a list of entries in preorder
	 */
	public ArrayList<String> list() {
		return root.list();
	}

	/**
	 * Find the nearest point to center. Ties are broken
	 * lexicographically.
	 */
	public LPoint nearestNeighbor(Point2D center) {
		LPoint best = null;
		LPoint nn = root.nearNeigh(center, rootCell, best, null);
		return nn;
	}

	/**
	 * Same as nearestNeighbor, but we return the nodes visited.
	 */
	public ArrayList<LPoint> nearestNeighborVisit(Point2D center) {
		LPoint best = null;
		ArrayList<LPoint> visited = new ArrayList<LPoint>();
		root.nearNeigh(center, rootCell, best, visited);
		Collections.sort(visited, new ByXThenY());
		return visited;
	}

	/**
	 * Generate an iterator.
	 */
	public LPointIterator iterator() {
		return new LPointIterator();
	}

	public LPoint closestCenter(LPoint point) {
		LPoint closeCenter =  root.closestCenter(point);
		return closeCenter;

	}

//	public <LPoint extends LabeledPoint2D> ArrayList<LPoint> getSites(LabeledPoint2D center, Rectangle2D cell, ArrayList<LPoint> list) {
//		ArrayList<LPoint> sites = new ArrayList<>();
//		return root.getSites(center, cell, sites);
//	}


	// -----------------------------------------------------------------
	// Tree Utilities for sorting and splitting.
	// -----------------------------------------------------------------

	/**
	 * Create a new subtree from a list of points. This does most
	 * of the work of the sliding midpoint method. If there is at
	 * most point point, we create a new external node. Otherwise,
	 * we determine the longest dimension of the current cell
	 * (breaking ties in favor of x) and make this the cutting
	 * dimension. We then compute the cell's midpoint and make this
	 * the cutting value. We sort the points by the cutting
	 * dimension and partition them about the cutting
	 * value. Points that lie on the cutting value are placed in
	 * the right subtree. If the resulting partition is trivial, we
	 * set the cutting value to the coordinate (largest or smallest).
	 */
	Node bulkCreate(List<LPoint> pts, Rectangle2D cell) throws Exception {
		Node u = bulkCreateHelper(pts, cell); // build a new tree

		// root contender?
		for(LPoint curr: root.contenders){
			u.addCenter(curr, cell);
		}
		return u;
	}

	Node bulkCreateHelper(List<LPoint> pts, Rectangle2D cell) {
		if (pts.size() == 0) { // no points
			return new ExternalNode(null);
		} else if (pts.size() == 1) { // one point
			return new ExternalNode(pts.get(0));
		} else { // two or more points
			int cutDim = getWideSide(cell); // get cutting dimension

			sortPts(pts, cutDim); // sort along cutting dimension
			LPoint leftmost = pts.get(0); // leftmost/rightmost points
			LPoint rightmost = pts.get(pts.size()-1);
			if (leftmost.get(cutDim) == rightmost.get(cutDim)) {
				cutDim = 1 - cutDim; // reverse splitting direction
			}
			double cutVal = (cell.low.get(cutDim) + cell.high.get(cutDim))/2; // default cutting value
			if (rightmost.get(cutDim) < cutVal) { // all points on left?
				cutVal = rightmost.get(cutDim);
			} else if (leftmost.get(cutDim) > cutVal) { // all points on right?
				cutVal = leftmost.get(cutDim);
			}
			int splitIndex = 0; // where to split
			while (splitIndex < pts.size() && pts.get(splitIndex).get(cutDim) < cutVal) {
				splitIndex++;
			}
			// left/right cells
			Rectangle2D leftCell = cell.leftPart(cutDim, cutVal);
			Rectangle2D rightCell = cell.rightPart(cutDim, cutVal);
			// the new node
			InternalNode u = new InternalNode(cutDim, cutVal);
			u.left = bulkCreateHelper(pts.subList(0, splitIndex), leftCell);
			u.right = bulkCreateHelper(pts.subList(splitIndex,  pts.size()), rightCell);
			u.left.parent = u.right.parent = u;
			u.size = u.left.getSize() + u.right.getSize();
			return u;
		}
	}

	/**
	 * Get the wider side of a rectangle.
	 */
	private static int getWideSide(Rectangle2D cell) {
		return (cell.getWidth(0) >= cell.getWidth(1) ? 0 : 1);
	}



	/**
	 * Return a reference to the closer of two points. If the two
	 * distances are the same, it returns the point having the smaller
	 * (x,y) coordinates in lexicographical order.
	 */
	private LPoint closerOf(Point2D q, LPoint p1, LPoint p2) {
		double dist1 = (p1 == null ? Double.POSITIVE_INFINITY : q.distanceSq(p1.getPoint2D()));
		double dist2 = (p2 == null ? Double.POSITIVE_INFINITY : q.distanceSq(p2.getPoint2D()));
		if      (p1 == null || dist1 > dist2) return p2;
		else if (p2 == null || dist2 > dist1) return p1;
		else { // same distances and both are non-null
			if      (p1.getX() < p2.getX()) return p1;
			else if (p2.getX() < p1.getX()) return p2;
			else { // ...and same x-coordinates
				if (p1.getY() <= p2.getY()) return p1;
				else return p2;
			}
		}
	}


	/**
	 * Sort points according to given dimension.
	 */
	private void sortPts(List<LPoint> pts, int d) {
		if (d == 0) { // sort by the cutting dimension
			Collections.sort(pts, new ByXThenY());
		} else {
			Collections.sort(pts, new ByYThenX());
		}
	}

	/**
	 * Comparators for sorting (by X or by Y). Both are lexicographical.
	 */
	private class ByXThenY implements Comparator<LPoint> { // lexicographic (x,y)
		public int compare(LPoint pt1, LPoint pt2) {
			double x1 = pt1.getX();
			double x2 = pt2.getX();
			if (x1 < x2)
				return -1;
			else if (x1 > x2)
				return +1;
			else {
				double y1 = pt1.getY();
				double y2 = pt2.getY();
				if (y1 < y2)
					return -1;
				else if (y1 > y2)
					return +1;
				else {
					return 0;
				}
			}
		}
	}

	private class ByYThenX implements Comparator<LPoint> { // lexicographic (y,x)
		public int compare(LPoint pt1, LPoint pt2) {
			double y1 = pt1.getY();
			double y2 = pt2.getY();
			if (y1 < y2)
				return -1;
			else if (y1 > y2)
				return +1;
			else {
				double x1 = pt1.getX();
				double x2 = pt2.getX();
				if (x1 < x2)
					return -1;
				else if (x1 > x2)
					return +1;
				else {
					return 0;
				}
			}
		}
	}

	class LabelComparator implements Comparator<LPoint> {
		public int compare(LPoint p1, LPoint p2) {
			return p1.getLabel().compareTo(p2.getLabel());
		}
	}
	
	// ----------------------------------------------------------------
	// Possible new additions for Programming Assignment 3
	// You may modify/add/remove these as you see fit, since they will 
	// just be used internally by your code.
	// ----------------------------------------------------------------
	 ArrayList<LPoint> addCenter(LPoint center, Rectangle2D cell) throws Exception {
		return root.addCenter(center, rootCell);
	}
	 ArrayList<String> listWithCenters() { return root.listWithCenters(); }

	public ArrayList<AssignedPair> listAssignments() {
		ArrayList<AssignedPair> pairs = root.listAssignments();
		Collections.sort(pairs, new Comparator<AssignedPair>() {
			@Override
			public int compare(AssignedPair o1, AssignedPair o2) {
				if (o1.distanceSq < o2.distanceSq) {
					return -1;
				} else if (o1.distanceSq > o2.distanceSq) {
					return 1;
				} else {
					double x1 = o1.site.getX();
					double x2 = o2.site.getX();
					if (x1 < x2)
						return -1;
					else if (x1 > x2)
						return +1;
					else {
						double y1 = o1.site.getY();
						double y2 = o2.site.getY();
						if (y1 < y2)
							return -1;
						else if (y1 > y2)
							return +1;
						else {
							return 0;
						}
					}
				}
			}
		});
		return pairs;
	}

}
