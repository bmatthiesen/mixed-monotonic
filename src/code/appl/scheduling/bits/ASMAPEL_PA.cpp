/* Copyright (C) 2018-2019 Bho Matthiesen, Christoph Hellings
 * 
 * This program is used in the article:
 *
 * Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang
 * Utschick, "Mixed Monotonic Programming for Fast Global Optimization,"
 * submitted to IEEE  Transactions on Signal Processing.
 * 
 * 
 * License:
 * This program is licensed under the GPLv2 license. If you in any way use this
 * code for research that results in publications, please cite our original
 * article listed above.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. */


#define CALL(F, X) ((F)(shift(X), user_data))

template <size_t Dim>
void
PA<Dim>::optimize(bool startPointPresent, typename PA<Dim>::vtype startPoint)
{
	VertexSet vertexSet;
	vtype b;

	std::chrono::time_point<clock> tic(clock::now());

	// shift
	if (doshift)
	{
		if (std::all_of(xshift.begin(), xshift.end(), [] (auto e) { return e == 0.0; }))
			xshift = ub;

		if (startPointPresent)
		{
			for (size_t i = 0; i < Dim; i++)
				startPoint[i] = startPoint[i] + xshift[i];
		}

		for (size_t i = 0; i < Dim; i++)
			b[i] = ub[i] + xshift[i];
	}
	else
		b = ub;

	// init xopt / optval
	if (startPointPresent)
	{
		if (CALL(inG, startPoint) && CALL(inH, startPoint))
		{
			xopt = startPoint;
			optval = CALL(obj, xopt);
		}
		else
			std::cerr << "Infeasible startPoint. Starting from scratch." << std::endl;

	}
	else
	{
		xopt = b;
		optval = -std::numeric_limits<double>::infinity();
	}

	// init other
	setStatus(UNSOLVED);
	iter = 0;
	lastUpdate = 0;

	vertexSet.emplace_back(b, CALL(obj, b));
	max_queue_size = vertexSet.size();

	while (true)
	{
		// find new candidate
		typename VertexSet::const_iterator zk = std::max_element(vertexSet.begin(), vertexSet.end(),
				[](const Vertex& a, const Vertex& b) { return a.obj < b.obj; }
			);

		// output
		if (output && iter != 0 && iter % output_every == 0)
			std::printf("%8llu %8lu %11g %11g (D = %g) | Peak RSS: %zu\n", iter, vertexSet.size(), optval, zk->obj, zk->obj - optval, getPeakRSS());

		max_queue_size = std::max(max_queue_size, vertexSet.size());

		// increase counter & check MaxIter
		if (++iter > MaxIter) // this does not break if MaxIter = (unsigned long long) -1 
		{
			setStatus(MAXITER);
			break;
		}

		// we're done if candidate point is feasible
		if (CALL(inG, zk->v))
		{
			xopt = zk->v;
			optval = zk->obj;
			lastUpdate = iter;
			setStatus(OPTIMAL);
			break;
		}

		Vertex zkV = *zk;
		typename VertexSet::iterator Tacc = partition(vertexSet.begin(), vertexSet.end(),
				[&](const Vertex& v)->bool // return true if v not in Zn
				{
					return ((1+tol)*v.obj < zkV.obj);
				}
			);

		VertexSet Zn;
		Zn.reserve(std::distance(Tacc, vertexSet.end()));

		Zn.insert(
				Zn.end(),
				std::move_iterator<typename VertexSet::iterator>(Tacc),
				std::move_iterator<typename VertexSet::iterator>(vertexSet.end())
			);

		vertexSet.erase(Tacc, vertexSet.end());
		
		// projection
		auto ret = projection(zkV.v, user_data);
		Vertex newPoint(std::move(ret.ret));
		Vertex newFeasiblePoint(std::move(ret.cand));

		// is point feasible (i.e. also in H)
		if (CALL(inH, newFeasiblePoint.v))
		{
			newFeasiblePoint.obj = CALL(obj, newFeasiblePoint.v);
			if (newFeasiblePoint.obj >= optval)
			{
				xopt = newFeasiblePoint.v;
				optval = newFeasiblePoint.obj;
				lastUpdate = iter;
			}
		}

		for (const auto& e : Zn)
		{
			// projection
			auto ret = projection(e.v, user_data);
			Vertex tmpPoint(std::move(ret.ret));
			Vertex tmpFeasiblePoint(std::move(ret.cand));

			// is point feasible (i.e. also in H)
			if (CALL(inH, tmpFeasiblePoint.v))
			{
				tmpFeasiblePoint.obj = CALL(obj, tmpFeasiblePoint.v);
				
				if (tmpFeasiblePoint.obj >= optval)
				{
					xopt = tmpFeasiblePoint.v;
					optval = tmpFeasiblePoint.obj;
					lastUpdate = iter;
				}
				
				if ((e.obj-tmpFeasiblePoint.obj) < (zkV.obj-newFeasiblePoint.obj))
				{
					newPoint = tmpPoint;
					newFeasiblePoint = tmpFeasiblePoint;
					zkV = e;
				}

			}
		}

		vertexSet.push_back(std::move(zkV));

		// compute T*
		typename VertexSet::iterator Tsit = partition(vertexSet.begin(), vertexSet.end(),
				[&](const Vertex& v)->bool // return true if v not in Tstar
				{
					for (size_t i = 0; i < Dim; i++)
					{
						if (v.v[i] <= newPoint.v[i])
						{
							return true;
						}
					}

					return false;
				}
			);

		VertexSet Tstar;
		Tstar.reserve(std::distance(Tsit, vertexSet.end()));

		Tstar.insert(
				Tstar.end(),
				std::move_iterator<typename VertexSet::iterator>(Tsit),
				std::move_iterator<typename VertexSet::iterator>(vertexSet.end())
			);

		vertexSet.erase(Tsit, vertexSet.end()); // TODO this should be done later

		// delete dominated
		vertexSet.erase(
				std::remove_if(
					vertexSet.begin(),
					vertexSet.end(),
					[this](const Vertex& v) { return v.obj <= optval; }
				),
				vertexSet.end()
			);

		// create new Vertices
		VertexSet newVertices;
		newVertices.reserve(Tstar.size() * Dim);
		for (typename VertexSet::const_iterator e = Tstar.begin(); e != Tstar.end(); e++)
		{
			for (size_t i = 0; i < Dim; i++)
			{
				Vertex v(e->v);
				v.v[i] = newPoint.v[i];

				if (CALL(inH, v.v))
				{
					v.obj = CALL(obj, v.v);

					if (v.obj > optval + epsilon)
						newVertices.push_back(std::move(v));
				}
			}
		}

		// remove improper new vertices
		newVertices.erase(
				std::remove_if(
					newVertices.begin(),
					newVertices.end(),
					[newVertices,this](const Vertex& v) // TODO is there a smarter way instead of passing newVertices as value???
					{
						for(auto const &e : newVertices)
						{
							bool le(true), neq(false);

							for (size_t i = 0; i < Dim; i++)
							{
								if (v.v[i] > e.v[i])
								{
									le = false;
									break;
								}
							}

							for (size_t i = 0; i < Dim; i++)
							{	
								if (v.v[i] != e.v[i])
								{
									neq = true;
									break;
								}
							}

							if (le && neq)
								return true;
						}

						return false;
					}
				),
				newVertices.end()
			);

		// add newVertices to vertexSet
		if (newVertices.size() > 0)
		{
			vertexSet.reserve(vertexSet.size() + newVertices.size());
			vertexSet.insert(
					vertexSet.end(),
					std::move_iterator<typename VertexSet::iterator>(newVertices.begin()),
					std::move_iterator<typename VertexSet::iterator>(newVertices.end())
				);
			newVertices.clear();
		}

		if (vertexSet.size() == 0)
		{
			setStatus(OPTIMAL);
			break;
		}
	}

	// undo shift
	if (doshift)
		xopt = shift(xopt);

	runtime = std::chrono::duration_cast<std::chrono::nanoseconds> (clock::now() - tic).count() / 1e9;

	// final status output
	if (output)
	{
		std::printf("%8llu %8lu %11g %11g (D = %g) | Peak RSS: %zu\n", iter, vertexSet.size(), optval, optval, 0.0, getPeakRSS());
		std::cout << std::endl;
		printResult();
	}
}

template <size_t Dim>
void
PA<Dim>::printResult()
{
	std::cout << "Status: " << statusStr << std::endl;
	std::cout << "Optval: " << optval << std::endl;

	std::cout << "X*: [";
	std::for_each(xopt.begin(), xopt.end(), [] (const basetype &a) { std::cout << " " << a; });
	std::cout << " ]" << std::endl;

	std::cout << "Precision: epsilon = " << epsilon << std::endl;

	std::cout << "Iter: " << iter << std::endl;
	std::cout << "Solution found in iter: " << lastUpdate << std::endl;

	std::cout << "Runtime: " << runtime << " sec" << std::endl;
	std::cout << "Peak RSS: " << getPeakRSS() << std::endl;
}

template <size_t Dim>
void
PA<Dim>::setShift(basetype e)
{
	vtype v;

	for (auto &el : v)
		el = e;

	setShift(v);
}
