"Some lazy decorators for methods of classes"
import numpy as np, os, functools

class ClassWraps(type):	# metaclass
	"""kinda like functools.wraps decorator, but a metaclass and specifically 
	for descriptor-decorator classes;
	retain the documentation of the decorated methods inside a class"""
	def __new__(cls, clsname, parents, attrs):
		def copy_attr(attr):	# closure
			return property(lambda self: getattr(self._method, attr))
		for attr in ['__name__', '__doc__', '__module__', '__annotations__']:	# `__qualname__` won't work like this...
			attrs[attr] = copy_attr(attr)
			# damn, `__name__` attribute is especially tricky here... 
			# but heck it works this way, so no whining!
			# (I managed to make it work through a metaclass forbidden magic)
			# this breaks documentation of the classes useing this metaclass,
			# but I don't care!
		return super().__new__(cls, clsname, parents, attrs)

class MethodDecorator(metaclass=ClassWraps):	# base class
	"""base class for descriptors which decorate methods; retains documentation of the methods;
	the original method is supposed to be accessed through the `_method` attribute;
	the internal `_obj` attribute is updated each time the descriptor is invoked through `__get__`;
	this docstring is lost in oblivion because of the metaclass"""
	def __init__(self, method):
		self._method = method
	def __get__(self, obj, objtype=None):
		self._obj = obj
		return self

class lazy_property(MethodDecorator):	# descriptor-decorator
	"""very lazy property decorator, just like myself;
	computes the value only once, and only when it is explicitly asked to;
	this docstring is lost, `help` is helpless this time"""
	def __get__(self, obj, objtype=None):
		if obj is None:
			return self
		value = self._method(obj)
		setattr(obj, self._method.__name__, value)
		return value

def get_dump_filename(filename, method_name):
	return filename + '.' + method_name + '.npy'

class _LazyArray:
	def __init__(self, instance, method, length, cached, dumped):
		self.instance = instance
		self.method = method
		self.length = length
		self.dumped = dumped
		self.func = lambda i: method(instance, i)
		if cached:
			self.func = functools.lru_cache(length)(self.func)
	@lazy_property
	def arr(self):
		arr = np.array([self.func(i) for i in range(self.length)])
		if self.dumped and len(arr) > 1:
			np.save(get_dump_filename(self.instance.filename, self.method.__name__), arr)
		setattr(self.instance, self.method.__name__, arr)	# quite dirty afterwards
		return arr
	def __getitem__(self, i):
		return self.func(i) if isinstance(i, int) else self.arr[i]
	def __iter__(self):
		return iter(self.arr)
	def __len__(self):
		return self.length
	def __getattr__(self, name):
		return getattr(self.arr, name)

def lazy_array(length=None, cached=True, dumped=False):	# closure of a descriptor-decorator
	class LazyArrayProperty(MethodDecorator):
		def __get__(self, obj, objtype=None):
			if obj is None:
				return self
			if dumped:
				dump_filename = get_dump_filename(obj.filename, self._method.__name__)
			# Note: in the following, checking modification time is not very reliable (when copied and not preserved)
			if dumped and os.path.isfile(dump_filename) and os.path.getmtime(dump_filename) > os.path.getmtime(obj.filename):
				# subject to race conditions, but I don't care
				value = np.load(dump_filename)
			else:
				_length = len(obj) if callable(length) or length is None else length
				value = _LazyArray(obj, self._method, _length, cached=cached, dumped=dumped)
			setattr(obj, self._method.__name__, value)
			return value
	if callable(length):
		return LazyArrayProperty(length)
	return LazyArrayProperty

# def lazy_array(length=None, cached=True, dumped=False):	# closure of a descriptor-decorator
# 	class LazyArrayProperty(MethodDecorator):
# 		def __init__(self, method):
# 			super().__init__(method)
# 			self._dumped = dumped
# 		def __get__(self, obj, objtype=None):
# 			self._length = len(self._obj) if callable(length) or length is None else length
# 			self._func = lambda i: self._method(self._obj, i)
# 			if cached:
# 				self._func = functools.lru_cache(self._length)(self._func)
# 			return super().__get__(obj, objtype)
# 		@lazy_property
# 		def _arr(self):
# 			dump_filename = get_dump_filename(self._obj.filename, self._method.__name__)
# 			if (self._dumped and os.path.isfile(dump_filename) and
# 				os.path.getmtime(dump_filename) > os.path.getmtime(self._obj.filename)):
# 				arr = np.load(dump_filename)
# 			else:
# 				arr = np.array([self._func(i) for i in range(self._length)])
# 				if self._dumped and len(arr) > 1:
# 					np.save(get_dump_filename(self._obj.filename, self._method.__name__), arr)
# 			return arr
# 		def __getitem__(self, i):
# 			return self._func(i) if isinstance(i, int) else self._arr[i]
# 		def __iter__(self):
# 			return iter(self._arr)
# 		def __len__(self):
# 			return self._length
# 		def __getattr__(self, name):
# 			return getattr(self._arr, name)
# 	if callable(length):
# 		return LazyArrayProperty(length)
# 	return LazyArrayProperty

