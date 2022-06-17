(function() {var implementors = {};
implementors["matrixmultiply"] = [{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"enum\" href=\"matrixmultiply/enum.CGemmOption.html\" title=\"enum matrixmultiply::CGemmOption\">CGemmOption</a>","synthetic":true,"types":["matrixmultiply::gemm::CGemmOption"]}];
implementors["ndarray"] = [{"text":"impl&lt;A&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.OwnedRepr.html\" title=\"struct ndarray::OwnedRepr\">OwnedRepr</a>&lt;A&gt;","synthetic":true,"types":["ndarray::data_repr::OwnedRepr"]},{"text":"impl&lt;'a, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.Axes.html\" title=\"struct ndarray::iter::Axes\">Axes</a>&lt;'a, D&gt;","synthetic":true,"types":["ndarray::dimension::axes::Axes"]},{"text":"impl&lt;D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.Indices.html\" title=\"struct ndarray::iter::Indices\">Indices</a>&lt;D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::indexes::Indices"]},{"text":"impl&lt;D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.IndicesIter.html\" title=\"struct ndarray::iter::IndicesIter\">IndicesIter</a>&lt;D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::indexes::IndicesIter"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.AxisChunksIter.html\" title=\"struct ndarray::iter::AxisChunksIter\">AxisChunksIter</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::AxisChunksIter"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.AxisChunksIterMut.html\" title=\"struct ndarray::iter::AxisChunksIterMut\">AxisChunksIterMut</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::AxisChunksIterMut"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.AxisIter.html\" title=\"struct ndarray::iter::AxisIter\">AxisIter</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::AxisIter"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.AxisIterMut.html\" title=\"struct ndarray::iter::AxisIterMut\">AxisIterMut</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::AxisIterMut"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.ExactChunks.html\" title=\"struct ndarray::iter::ExactChunks\">ExactChunks</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::chunks::ExactChunks"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.ExactChunksIter.html\" title=\"struct ndarray::iter::ExactChunksIter\">ExactChunksIter</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::chunks::ExactChunksIter"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.ExactChunksIterMut.html\" title=\"struct ndarray::iter::ExactChunksIterMut\">ExactChunksIterMut</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::chunks::ExactChunksIterMut"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.ExactChunksMut.html\" title=\"struct ndarray::iter::ExactChunksMut\">ExactChunksMut</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::chunks::ExactChunksMut"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.IndexedIter.html\" title=\"struct ndarray::iter::IndexedIter\">IndexedIter</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::IndexedIter"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.IndexedIterMut.html\" title=\"struct ndarray::iter::IndexedIterMut\">IndexedIterMut</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::IndexedIterMut"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.Iter.html\" title=\"struct ndarray::iter::Iter\">Iter</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::Iter"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.IterMut.html\" title=\"struct ndarray::iter::IterMut\">IterMut</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::IterMut"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.Lanes.html\" title=\"struct ndarray::iter::Lanes\">Lanes</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::lanes::Lanes"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.LanesIter.html\" title=\"struct ndarray::iter::LanesIter\">LanesIter</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::LanesIter"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.LanesIterMut.html\" title=\"struct ndarray::iter::LanesIterMut\">LanesIterMut</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::LanesIterMut"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.LanesMut.html\" title=\"struct ndarray::iter::LanesMut\">LanesMut</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::lanes::LanesMut"]},{"text":"impl&lt;'a, A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/iter/struct.Windows.html\" title=\"struct ndarray::iter::Windows\">Windows</a>&lt;'a, A, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::iterators::windows::Windows"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.ShapeError.html\" title=\"struct ndarray::ShapeError\">ShapeError</a>","synthetic":true,"types":["ndarray::error::ShapeError"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"enum\" href=\"ndarray/enum.ErrorKind.html\" title=\"enum ndarray::ErrorKind\">ErrorKind</a>","synthetic":true,"types":["ndarray::error::ErrorKind"]},{"text":"impl&lt;T&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.MathCell.html\" title=\"struct ndarray::MathCell\">MathCell</a>&lt;T&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;T: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::math_cell::MathCell"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"enum\" href=\"ndarray/enum.Order.html\" title=\"enum ndarray::Order\">Order</a>","synthetic":true,"types":["ndarray::order::Order"]},{"text":"impl&lt;D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.Shape.html\" title=\"struct ndarray::Shape\">Shape</a>&lt;D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::shape_builder::Shape"]},{"text":"impl&lt;D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.StrideShape.html\" title=\"struct ndarray::StrideShape\">StrideShape</a>&lt;D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::shape_builder::StrideShape"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.Slice.html\" title=\"struct ndarray::Slice\">Slice</a>","synthetic":true,"types":["ndarray::slice::Slice"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.NewAxis.html\" title=\"struct ndarray::NewAxis\">NewAxis</a>","synthetic":true,"types":["ndarray::slice::NewAxis"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"enum\" href=\"ndarray/enum.SliceInfoElem.html\" title=\"enum ndarray::SliceInfoElem\">SliceInfoElem</a>","synthetic":true,"types":["ndarray::slice::SliceInfoElem"]},{"text":"impl&lt;T, Din, Dout&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.SliceInfo.html\" title=\"struct ndarray::SliceInfo\">SliceInfo</a>&lt;T, Din, Dout&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;Din: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;Dout: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;T: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::slice::SliceInfo"]},{"text":"impl&lt;Parts, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.Zip.html\" title=\"struct ndarray::Zip\">Zip</a>&lt;Parts, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;Parts: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::zip::Zip"]},{"text":"impl&lt;T&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"enum\" href=\"ndarray/enum.FoldWhile.html\" title=\"enum ndarray::FoldWhile\">FoldWhile</a>&lt;T&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;T: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::zip::FoldWhile"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.AxisDescription.html\" title=\"struct ndarray::AxisDescription\">AxisDescription</a>","synthetic":true,"types":["ndarray::dimension::axes::AxisDescription"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.Axis.html\" title=\"struct ndarray::Axis\">Axis</a>","synthetic":true,"types":["ndarray::dimension::axis::Axis"]},{"text":"impl&lt;I:&nbsp;?<a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Sized.html\" title=\"trait core::marker::Sized\">Sized</a>&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.Dim.html\" title=\"struct ndarray::Dim\">Dim</a>&lt;I&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;I: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::dimension::dim::Dim"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.IxDynImpl.html\" title=\"struct ndarray::IxDynImpl\">IxDynImpl</a>","synthetic":true,"types":["ndarray::dimension::dynindeximpl::IxDynImpl"]},{"text":"impl&lt;S, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.ArrayBase.html\" title=\"struct ndarray::ArrayBase\">ArrayBase</a>&lt;S, D&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;D: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,<br>&nbsp;&nbsp;&nbsp;&nbsp;S: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::ArrayBase"]},{"text":"impl&lt;A&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.OwnedArcRepr.html\" title=\"struct ndarray::OwnedArcRepr\">OwnedArcRepr</a>&lt;A&gt;","synthetic":true,"types":["ndarray::OwnedArcRepr"]},{"text":"impl&lt;A&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.RawViewRepr.html\" title=\"struct ndarray::RawViewRepr\">RawViewRepr</a>&lt;A&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;A: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::RawViewRepr"]},{"text":"impl&lt;A&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"ndarray/struct.ViewRepr.html\" title=\"struct ndarray::ViewRepr\">ViewRepr</a>&lt;A&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;A: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["ndarray::ViewRepr"]},{"text":"impl&lt;'a, A&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"enum\" href=\"ndarray/enum.CowRepr.html\" title=\"enum ndarray::CowRepr\">CowRepr</a>&lt;'a, A&gt;","synthetic":true,"types":["ndarray::CowRepr"]}];
implementors["num_complex"] = [{"text":"impl&lt;T&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"num_complex/struct.Complex.html\" title=\"struct num_complex::Complex\">Complex</a>&lt;T&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;T: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["num_complex::Complex"]},{"text":"impl&lt;E&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"num_complex/struct.ParseComplexError.html\" title=\"struct num_complex::ParseComplexError\">ParseComplexError</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;E: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["num_complex::ParseComplexError"]}];
implementors["num_integer"] = [{"text":"impl&lt;A&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"num_integer/struct.ExtendedGcd.html\" title=\"struct num_integer::ExtendedGcd\">ExtendedGcd</a>&lt;A&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;A: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["num_integer::ExtendedGcd"]},{"text":"impl&lt;T&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"num_integer/struct.IterBinomial.html\" title=\"struct num_integer::IterBinomial\">IterBinomial</a>&lt;T&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;T: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a>,&nbsp;</span>","synthetic":true,"types":["num_integer::IterBinomial"]}];
implementors["num_traits"] = [{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"enum\" href=\"num_traits/enum.FloatErrorKind.html\" title=\"enum num_traits::FloatErrorKind\">FloatErrorKind</a>","synthetic":true,"types":["num_traits::FloatErrorKind"]},{"text":"impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.61.0/core/marker/trait.Unpin.html\" title=\"trait core::marker::Unpin\">Unpin</a> for <a class=\"struct\" href=\"num_traits/struct.ParseFloatError.html\" title=\"struct num_traits::ParseFloatError\">ParseFloatError</a>","synthetic":true,"types":["num_traits::ParseFloatError"]}];
if (window.register_implementors) {window.register_implementors(implementors);} else {window.pending_implementors = implementors;}})()