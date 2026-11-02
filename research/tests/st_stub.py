"""Streamlit 없이 app.py 를 끝까지 실행해 보는 스모크 테스트용 가짜 streamlit 모듈.

위젯은 기본값을 돌려주고, OVERRIDES[label] 로 특정 위젯 값을 지정할 수 있다.
그림(st.pyplot)은 FIGS 에 모아 두었다가 저장할 수 있다.
"""
import sys
import types

OVERRIDES = {}
FIGS = []
CALLS = []


class _Stop(Exception):
    pass


class _Ctx:
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def __getattr__(self, name):
        return getattr(_api, name)


def _val(label, default):
    CALLS.append(label)
    return OVERRIDES.get(label, default)


class _API:
    session_state = {}

    def __getattr__(self, name):  # 표시 전용 함수 (title, markdown, info, ...)
        def f(*a, **k):
            CALLS.append(name)
            return None
        return f

    def set_page_config(self, **k): pass

    def cache_data(self, *a, **k):
        if a and callable(a[0]):
            return a[0]
        return lambda fn: fn

    def cache_resource(self, *a, **k):
        return self.cache_data(*a, **k)

    def rerun(self):
        raise _Stop()

    def stop(self):
        raise _Stop()

    def columns(self, spec, **k):
        n = spec if isinstance(spec, int) else len(spec)
        return [_Ctx() for _ in range(n)]

    def tabs(self, labels):
        return [_Ctx() for _ in labels]

    def spinner(self, *a, **k):
        return _Ctx()

    def expander(self, *a, **k):
        return _Ctx()

    def container(self, *a, **k):
        return _Ctx()

    def radio(self, label, options, index=0, **k):
        return _val(label, options[index])

    def selectbox(self, label, options, index=0, **k):
        return _val(label, list(options)[index])

    def multiselect(self, label, options, default=None, **k):
        return _val(label, default or [])

    def slider(self, label, min_value=None, max_value=None, value=None, step=None, **k):
        return _val(label, value)

    def select_slider(self, label, options, value=None, **k):
        return _val(label, value)

    def number_input(self, label, min_value=None, max_value=None, value=None, step=None, **k):
        return _val(label, value)

    def text_input(self, label, value="", **k):
        return _val(label, value)

    def file_uploader(self, label, **k):
        return _val(label, None)

    def checkbox(self, label, value=False, **k):
        return _val(label, value)

    def button(self, label, **k):
        return _val(label, False)

    def pyplot(self, fig, **k):
        FIGS.append(fig)


_api = _API()
mod = types.ModuleType("streamlit")
for _n in dir(_API):
    if not _n.startswith("__"):
        setattr(mod, _n, getattr(_api, _n))
mod.__getattr__ = lambda name: getattr(_api, name)
mod.sidebar = _Ctx()
mod.session_state = _api.session_state
mod.Stop = _Stop


comp = types.ModuleType("streamlit.components")
comp_v1 = types.ModuleType("streamlit.components.v1")
comp_v1.html = lambda *a, **k: CALLS.append("components.html")
comp.v1 = comp_v1
mod.components = comp


def install():
    sys.modules["streamlit"] = mod
    sys.modules["streamlit.components"] = comp
    sys.modules["streamlit.components.v1"] = comp_v1
    return mod
